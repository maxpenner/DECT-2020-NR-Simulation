classdef dect_rx < handle
    
    properties
        verbose;        % show data during execution: 0 false, 1 only text, 2 text + plots
        mac_meta;       % data received from MAC layer
        phy_4_5;        % data from chapter 4 and 5

        packet_data;    % intermediate results during packet decoding
        
        tx_handle;      % handle to tx for debugging
        ch_handle;      % handle to rf channel for debugging

        STF_templates;  % for synchronization
        harq_buf_40;    % these are the cbsbuffers used by matlab for harq operation, PCC 40 bits
        harq_buf_80;    % these are the cbsbuffers used by matlab for harq operation, PCC 80 bits
        harq_buf;       % these are the cbsbuffers used by matlab for harq operation, PDC
        
        wiener;         % structure with data about wiener filter for channel estimation
    end
    
    methods
        function obj = dect_rx(verbose_arg, mac_meta_arg)
            obj.verbose = verbose_arg;
            obj.mac_meta = mac_meta_arg;
            obj.phy_4_5 = lib_util.run_chapter_4_5(verbose_arg, mac_meta_arg);

            obj.packet_data = [];

            obj.tx_handle = [];
            obj.ch_handle = [];

            % we only load the templates if they are need at any stage within the receiver
            if mac_meta_arg.sto_config.use_sto_sync == true || mac_meta_arg.sto_config.fine.fractional_cfo == true
                obj.STF_templates = lib_rx.sync_STF_template(mac_meta_arg);
            end
            obj.harq_buf_40 = [];
            obj.harq_buf_80 = [];
            obj.harq_buf = [];
            
            obj.wiener = [];
        end
        
        % We pass on the the samples as they are received by the antennas and we try to extract the PCC and PDC bits.
        % It gets more complicated when using HARQ as calls to this function depend on each other.
        function [PCC_user_bits, PDC_user_bits] = demod_decode_packet(obj, samples_antenna_rx)

            % for the purpose of readability, extract all variables that are necessary at this stage

            verbose_            = obj.verbose;
            tx_handle_          = obj.tx_handle;
            ch_handle_          = obj.ch_handle;
            STF_templates_      = obj.STF_templates;
            harq_buf_40_        = obj.harq_buf_40;
            harq_buf_80_        = obj.harq_buf_80;
            harq_buf_           = obj.harq_buf;
            wiener_             = obj.wiener;
        
            mode_0_to_1         = obj.phy_4_5.tm_mode.mode_0_to_1;
            N_eff_TX            = obj.phy_4_5.tm_mode.N_eff_TX;

            modulation0         = obj.phy_4_5.mcs.modulation0;

            N_b_DFT             = obj.phy_4_5.numerology.N_b_DFT;            
            N_b_CP              = obj.phy_4_5.numerology.N_b_CP;

            N_TB_bits           = obj.phy_4_5.N_TB_bits;
            N_PACKET_symb       = obj.phy_4_5.N_PACKET_symb;
            k_b_OCC             = obj.phy_4_5.k_b_OCC;

            n_packet_samples    = obj.phy_4_5.n_packet_samples;

            u                   = obj.mac_meta.u;
            Z                   = obj.mac_meta.Z;
            network_id          = obj.mac_meta.network_id;
            PLCF_type           = obj.mac_meta.PLCF_type;
            rv                  = obj.mac_meta.rv;
            N_RX                = obj.mac_meta.N_RX;
            oversampling        = obj.mac_meta.oversampling;

            sto_config          = obj.mac_meta.sto_config;
            cfo_config          = obj.mac_meta.cfo_config;

            use_ch_estim_type   = obj.mac_meta.use_ch_estim_type;

            use_equalization    = obj.mac_meta.use_equalization;

            physical_resource_mapping_PCC_cell = obj.phy_4_5.physical_resource_mapping_PCC_cell;
            physical_resource_mapping_PDC_cell = obj.phy_4_5.physical_resource_mapping_PDC_cell;
            physical_resource_mapping_STF_cell = obj.phy_4_5.physical_resource_mapping_STF_cell;
            physical_resource_mapping_DRS_cell = obj.phy_4_5.physical_resource_mapping_DRS_cell;

            %% sync of STO and extraction of N_eff_TX
            if sto_config.use_sto_sync == true

                % The number of samples received is larger than the number of samples in a packet.
                % In this function, we try to synchronize the packet and extract the exact number of samples in the packet.
                % In a real receiver, the number of samples is unknown, we would first have to decode the PCC which lies in the first few OFDM symbols.
                [samples_antenna_rx_sto_cfo, STO_CFO_report] = lib_rx.sync_STO( verbose_,...
                                                                                u,...
                                                                                N_b_DFT,...
                                                                                samples_antenna_rx,...
                                                                                STF_templates_.time_domain,...
                                                                                n_packet_samples,...
                                                                                oversampling,...
                                                                                sto_config,...
                                                                                cfo_config);

                obj.packet_data.STO_report = STO_CFO_report;
            else
                % assume packet is synchronized and has the correct length
                samples_antenna_rx_sto_cfo = samples_antenna_rx;
            end
            
            %% fractional CFO post STO sync and pre FFT based on STF
            if cfo_config.use_cfo_fractional == true

                % for fractional CFO correction we don't need the shape of STF, only length and number of patterns (derived from u)
                n_STF_samples = numel(cell2mat(STF_templates_.time_domain(1)));

                [samples_antenna_rx_sto_cfo, CFO_report] = lib_rx.sync_CFO_fractional(  verbose_,...
                                                                                        samples_antenna_rx_sto_cfo,...
                                                                                        n_STF_samples,...
                                                                                        u,...
                                                                                        cfo_config);

                obj.packet_data.CFO_report = CFO_report;
            end

            %% OFDM demodulation a.k.a FFT
            % Switch back to frequency domain.
            % We use one version with subcarriers from oversampling removed and one with them still occupied.
            % The second version is necessary if we have a very large integer CFO.
            [antenna_streams_mapped_rev, antenna_streams_mapped_rev_with_os_carriers] = lib_6_generic_procedures.ofdm_signal_generation_Cyclic_prefix_insertion_rev(samples_antenna_rx_sto_cfo,...
                                                                                                                                                                    k_b_OCC,...
                                                                                                                                                                    N_PACKET_symb,...
                                                                                                                                                                    N_RX,...
                                                                                                                                                                    N_eff_TX,...
                                                                                                                                                                    N_b_DFT,...
                                                                                                                                                                    u,...
                                                                                                                                                                    N_b_CP,...
                                                                                                                                                                    oversampling);

            %% sync of integer CFO post FFT based on STF
            if cfo_config.use_cfo_integer == true

                % sanity check
                if cfo_config.use_cfo_fractional == false
                    error('Correcting integer CFO, but not fractional CFO makes only sense for debugging purposes.');
                end

                [antenna_streams_mapped_rev, CFO_report_integer] = lib_rx.sync_CFO_integer( antenna_streams_mapped_rev_with_os_carriers,...
                                                                                            STF_templates_.freq_domain,...
                                                                                            N_RX,...
                                                                                            N_eff_TX,...
                                                                                            N_b_DFT,...
                                                                                            oversampling,...
                                                                                            cfo_config);

                % add integer CFO to report
                obj.packet_data.CFO_report_integer = CFO_report_integer;
            else
                % do nothing
            end
                                                                                                                    
            %% PCC decoding and packet data extraction
            % Decode PCC and extract all relevant data about the packet.
            % We know N_eff_TX from STO sync, so we also know the DRS pattern.
            % With STF and DRS known, we can now estimate the channel for PCC, which lies at the beginning of each packet.
            % PCC is always transmitted with transmit diversity coding.
            % For now it is done down below together with the PDC.
            % TODO

            %% residual CFO post FFT based on STF and DRS
            if cfo_config.use_cfo_residual == true
                antenna_streams_mapped_rev = lib_rx.sync_CFO_residual(  antenna_streams_mapped_rev,...
                                                                        physical_resource_mapping_DRS_cell,...
                                                                        physical_resource_mapping_STF_cell,...
                                                                        N_RX,...
                                                                        N_eff_TX,...
                                                                        cfo_config);
            end
                                                                                                                    
            %% channel estimation

            % Now that we know the length of the packet from the PCC, we can determine a channel estimate.
            % Output is a cell(N_RX,1), each cell with a matrix of size N_b_DFT x N_PACKET_symb x N_TX.
            % For subcarriers which are unused the channel estimate can be NaN or +/- infinity.
            if strcmp(use_ch_estim_type,'perfect') == true

                % only for AWGN OR rayleigh/rician with doppler frequency of 0 and flat fading
                ch_estim = lib_rx.channel_estimation_perfect(antenna_streams_mapped_rev, N_RX, N_eff_TX, ch_handle_);

            elseif strcmp(use_ch_estim_type,'perfect SIMO') == true
                
                % this is a cheap trick: we take the time domain samples without noise, FFT and take the effective channel ceofficients as a channel estimate
                antenna_streams_mapped_rev_no_noise = lib_6_generic_procedures.ofdm_signal_generation_Cyclic_prefix_insertion_rev(  ch_handle_.samples_antenna_rx_no_noise,...
                                                                                                                                    k_b_OCC,...
                                                                                                                                    N_PACKET_symb,...
                                                                                                                                    N_RX,...
                                                                                                                                    N_eff_TX,...
                                                                                                                                    N_b_DFT,...
                                                                                                                                    u,...
                                                                                                                                    N_b_CP,...
                                                                                                                                    oversampling);                
                
                ch_estim = lib_rx.channel_estimation_perfect_SIMO(antenna_streams_mapped_rev_no_noise, N_RX, N_eff_TX, tx_handle_);  

            elseif strcmp(use_ch_estim_type,'least squares') == true

                % Zero forcing at the pilots and linear interpolation for subcarriers inbetween.
                % It works without channel knowledge, but has a fairly large error floor.
                % The error floor is comes from only considering the closest neighbours.
                ch_estim = lib_rx.channel_estimation_ls(antenna_streams_mapped_rev, physical_resource_mapping_STF_cell, physical_resource_mapping_DRS_cell, N_RX, N_eff_TX);

            elseif strcmp(use_ch_estim_type,'wiener') == true

                % real world channel estimation based in precalculated wiener filter coefficients
                ch_estim = lib_rx.channel_estimation_wiener(antenna_streams_mapped_rev, physical_resource_mapping_DRS_cell, wiener_, N_RX, N_eff_TX);

            else
                error('Unknown channel estimation type %s.', use_ch_estim_type);
            end
            
            %% equalization
            % With the known transmission mode and the channel estimate, we can now equalize the received data.
            % This step is also often referred to as 'symbol detection'.
            if use_equalization == true
                
                % SISO
                if N_eff_TX == 1 && N_RX ==1

                    x_PCC_rev = lib_rx.equalization_SISO_zf(antenna_streams_mapped_rev, ch_estim, physical_resource_mapping_PCC_cell);
                    x_PDC_rev = lib_rx.equalization_SISO_zf(antenna_streams_mapped_rev, ch_estim, physical_resource_mapping_PDC_cell);

                % SIMO (Maximum Ratio Combining (MRC))
                elseif N_eff_TX == 1 && N_RX > 1

                    x_PCC_rev = lib_rx.equalization_SIMO_mrc(antenna_streams_mapped_rev, ch_estim, physical_resource_mapping_PCC_cell, N_RX);
                    x_PDC_rev = lib_rx.equalization_SIMO_mrc(antenna_streams_mapped_rev, ch_estim, physical_resource_mapping_PDC_cell, N_RX);

                % MISO (Alamouti) and MIMO (Alamouti + MRC)
                else
                    % Transmit diversity precoding
                    if ismember(mode_0_to_1, [1,5,10]) == true
                        x_PCC_rev = lib_rx.equalization_MISO_MIMO_alamouti_mrc(antenna_streams_mapped_rev, ch_estim, N_RX, N_eff_TX, physical_resource_mapping_PCC_cell);
                        x_PDC_rev = lib_rx.equalization_MISO_MIMO_alamouti_mrc(antenna_streams_mapped_rev, ch_estim, N_RX, N_eff_TX, physical_resource_mapping_PDC_cell);
                    end

                    % MIMO modes with more than one spatial stream
                    % TODO
                end
                
            % no equalization, a great debugging tool when using awgn channel
            else
                x_PCC_rev = lib_7_Transmission_modes.subcarrier_unmapping_PCC(antenna_streams_mapped_rev, physical_resource_mapping_PCC_cell);
                x_PDC_rev = lib_7_Transmission_modes.subcarrier_unmapping_PDC(antenna_streams_mapped_rev, physical_resource_mapping_PDC_cell);                
            end

            %% PCC and PDC decoding
            % decode PCC
            [PCC_user_bits, PCC_harq_buf_40_report,...
                            PCC_harq_buf_80_report,...
                            CL_report,...
                            BF_report,...
                            PCC_report,...
                            pcc_dec_dbg] =  lib_7_Transmission_modes.PCC_decoding(  x_PCC_rev, ...
                                                                                    harq_buf_40_,...
                                                                                    harq_buf_80_);
            
            % decode PDC
            [PDC_user_bits, PDC_harq_buf_report, pdc_dec_dbg] = lib_7_Transmission_modes.PDC_decoding(  x_PDC_rev,...
                                                                                                        N_TB_bits,...
                                                                                                        Z,...
                                                                                                        network_id,...
                                                                                                        PLCF_type,...
                                                                                                        rv,...
                                                                                                        modulation0,...
                                                                                                        harq_buf_);

            %% HARQ buffers
            % save harq buffer for next call
            if numel(PCC_harq_buf_40_report) > 0
                obj.harq_buf_40 = PCC_harq_buf_40_report;
            end
            if numel(PCC_harq_buf_80_report) > 0
                obj.harq_buf_80 = PCC_harq_buf_80_report;
            end            
            if numel(PDC_harq_buf_report) > 0
                obj.harq_buf = PDC_harq_buf_report;
            end

            %% save packet data for debugging
            obj.packet_data.CL_report   = CL_report;
            obj.packet_data.BF_report   = BF_report;
            obj.packet_data.PCC_report  = PCC_report;
            obj.packet_data.pcc_dec_dbg = pcc_dec_dbg;
            obj.packet_data.pdc_dec_dbg = pdc_dec_dbg;

            % debugging
            if verbose_ > 0
                disp('##### RX debugging ######');
                
                % compare complex symbols of PCC and PDF
                fprintf('Max error distance PCC: %.20f\n', max(abs(x_PCC_rev - tx_handle_.packet_data.x_PCC)));
                fprintf('Max error distance PDC: %.20f\n', max(abs(x_PDC_rev - tx_handle_.packet_data.x_PDC)));

                % measure the sinr in dB by comparing the complex PDC symbols
                tx_power_measurement = tx_handle_.packet_data.x_PDC;
                tx_power_measurement = sum(abs(tx_power_measurement).^2)/numel(tx_power_measurement);
                noise_power_measurement = x_PDC_rev - tx_handle_.packet_data.x_PDC;
                noise_power_measurement = sum(abs(noise_power_measurement).^2)/numel(noise_power_measurement);
                sinr = 10*log10(tx_power_measurement/noise_power_measurement);
                fprintf('Measured f domain SINR: %f dB\n', sinr);

                if verbose_ > 1
                    scatterplot(x_PDC_rev);
                end
            end
        end
    end
end


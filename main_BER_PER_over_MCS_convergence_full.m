clear all;
close all;

% This scripts can be used to test aspects which are usually not part of link-level simulations:
%
%      - resampling at TX (we assume USRP sample rates, master clock is a multiple of 30.72 MHz)
%      - ADC resolution at RX (after rf channel and noise)
%      - resampling at RX (includes narrow noise filter)
%      - u detection
%      - STO synchronization (based on STF)
%      - fractional CFO synchronization (based on STF, corrects into +/-2% range of subcarrier spacing)
%      - integer CFO synchronization (based on STF, can correct 60ppm at 6GHZ and 27kHz subcarrier spacing)
%      - b detection
%      - residual CFO and SCO synchronization (based on DRS, residual offset due to Carrier Frequency Offset (CFO) and Sample Clock Offset (SCO))
%      - Wiener Filter Channel Estimation (also part of Link-Level Simulations)

rng('shuffle');
%rng(1140598280);

%warning('off');

if exist('results', 'dir')
    lib_util.clear_directory('results');
else
    mkdir('results');
end

fprintf('Starting at %s\n', datestr(now,'HH:MM:SS'));

%profile on

% choose mcs to simulate and maximum number of harq retransmissions
mcs_index_vec = [1];
max_harq_retransmissions = 0;

% do we resamples to USRP rates?
use_resampling = true;
if use_resampling == true

    % the USRP samp rate should be at least this factor times the oversampled nominal DECT-2020 rate
    USRP_samp_rate_min_multiple = 1;
end

% quantization parameter so simulate ADC and low noise condition
use_ADC_scaling = true;
if use_ADC_scaling == true

    % force signal into ADC range after RF channel
    %ADC_range_scalingFactor = 0.7;     % severe clipping, but RX is still functional
    ADC_range_scalingFactor = 0.25;  	% avoid most clipping
    %ADC_range_scalingFactor = 0.01;    % very slow amplitude (make sure QPSK is used)
    
    % typical values are 8 (e.g. HackRF, RTL-SDR), 12, 14 and 16
    n_ADC_bits = 14;
end

% simulation range for link level simulation
snr_db_vec_global = 0:5:50;
snr_db_vec_global = repmat(snr_db_vec_global,numel(mcs_index_vec),1);

% convergence settings for per over snr packet simulation
n_packets_per_batch = 5e3;              % batch size for convergence comparison, this number should no be too low to avoid parfor-overhead
n_compare = 3;                          % how many batches do we compare?
compare_convergence_threshold = 0.01;   % maximum relative change from batch to batch
n_max_packets_snr = 0.50e6;             % we might never achieve convergence, in this case we have to limit simulation time
per_abort_threshold = 1e-5;             % we are not interested in PERs below this threshold
n_min_corruped_packets = 5;             % for very low PERs (<1e-4), a single packet can indicate convergence causing a bumby PER curve, therefore we need at least 5 erroneous packets

n_worker = 18;                                                  % how many workers? best case number this is the number of cpu-cores in your system (parfor)
n_packets_per_worker_call = ceil(n_packets_per_batch/n_worker); % how many packets does a single worker calculate per batch

% results container for STO synchronization
edges = [-10e6, -30.5 : 1 : 30.5, 10e6];
n_bins = numel(edges)-1;
STO_hist_coarse_global = zeros(n_bins, numel(snr_db_vec_global(1,:)));
STO_hist_fine_global = zeros(n_bins, numel(snr_db_vec_global(1,:)));

% result container for PCC
n_bits_PCC_sent_global = zeros(numel(mcs_index_vec), numel(snr_db_vec_global(1,:)));        % BER uncoded
n_bits_PCC_error_global = zeros(numel(mcs_index_vec), numel(snr_db_vec_global(1,:)));       % BER uncoded
n_packets_PCC_sent_global = zeros(numel(mcs_index_vec), numel(snr_db_vec_global(1,:)));     % PER
n_packets_PCC_error_global = zeros(numel(mcs_index_vec), numel(snr_db_vec_global(1,:)));    % PER

% result container for PDC
n_bits_PDC_sent_global = zeros(numel(mcs_index_vec), numel(snr_db_vec_global(1,:)));        % BER uncoded
n_bits_PDC_error_global = zeros(numel(mcs_index_vec), numel(snr_db_vec_global(1,:)));       % BER uncoded
n_packets_PDC_sent_global = zeros(numel(mcs_index_vec), numel(snr_db_vec_global(1,:)));     % PER
n_packets_PDC_error_global = zeros(numel(mcs_index_vec), numel(snr_db_vec_global(1,:)));    % PER

bps_global = zeros(numel(mcs_index_vec), 1);        % bits per symbol, PDC only
tbs_global = zeros(numel(mcs_index_vec), 1);        % transport block size, PDC only

cnt = 1;
for mcs_index = mcs_index_vec

    % these variables need to be set before creating tx and rx
    mac_meta_tx.u = 1;                                  % mu = 1, 2, 4 or 8
    mac_meta_tx.b = 1;                                  % beta = 1, 2, 4, 8, 12 or 16
    mac_meta_tx.PacketLengthType = 0;                   % 0 for subslots, 1 for slots
    mac_meta_tx.PacketLength = 2;                       % min is 1, max is 16 according to Table 6.2.1-2a in part 4
    mac_meta_tx.tm_mode_0_to_11 = 0;                 	% Table 7.2-1, mode determines wether transmission is closed loop or not, values range from 0 to 11
    mac_meta_tx.mcs_index = mcs_index;                  % Table A-1 in part 3, values range from 0 to 11
    mac_meta_tx.Z = 6144;                               % 5.3 -> so far only Z=6144 fully supported, 2048 only at TX, RX missing (Matlab has no option for Z=2048 in LTE toolbox)
    mac_meta_tx.oversampling = 2;                    	% By how much do we oversample our ofdm packet compared to critical sampling (insert zeros at specturm edges before IFFT)?
    mac_meta_tx.codebook_index = 0;                     % 6.3.4, any value other than 0 makes packet beamformed, throws error if out of bound (depends on tm_mode_0_to_11)
    mac_meta_tx.PLCF_type = 2;                          % Type 1 is 40 bits, Type 2 is 80 bits
    mac_meta_tx.rv = 0;                                 % HARQ version, values range from 0, 1, 2 to 3 (right HARQ retransmission order is 0 2 3 1)
    mac_meta_tx.network_id = de2bi(1e6,32,'left-msb');  % 7.6.6 must be given as a 32 bit vector with network_id(1) being the MSB, network_id must be known for scrambler on PHY 
    
    % temporary restrictions
    if mac_meta_tx.Z ~= 6144
        error('Z must be 6144.');
    end    
    
    % create tx
    verbose = 0;
    txx = dect_tx(verbose, mac_meta_tx);
    
    % determine optimal Hardware sample rate
    if use_resampling == true

        % find optimal USRP rate
        [USRP_samp_rate, L, M] = lib_rx.resampling_USRP_rate(txx.phy_4_5.numerology.B_u_b_DFT, mac_meta_tx.oversampling, USRP_samp_rate_min_multiple);

        % design resampling fir filters
        [fir_coef_tx, fir_coef_rx] = lib_rx.resampling_filter_design(txx.mac_meta.u, txx.mac_meta.b, txx.phy_4_5.numerology.N_b_DFT, txx.phy_4_5.numerology.n_guards_top, mac_meta_tx.oversampling, L, M);
    else

        % no resampling, required for RF channel
        L = 1;
        M = 1;
    end
    
    % additional rx configuration
    mac_meta_rx = mac_meta_tx;
    mac_meta_rx.N_RX = 1;

    % STO synchronization
    mac_meta_rx.sto_config = lib_rx.sync_STO_param(mac_meta_tx.u, mac_meta_tx.b, mac_meta_tx.oversampling);
    mac_meta_rx.sto_config.use_sto_sync = true;

    % CFO synchronization
    mac_meta_rx.cfo_config = lib_rx.sync_CFO_param(mac_meta_tx.u);
    mac_meta_rx.cfo_config.use_cfo_fractional = true;
    mac_meta_rx.cfo_config.use_cfo_integer = true;
    mac_meta_rx.cfo_config.use_cfo_residual = true;

    % channel estimation
    mac_meta_rx.use_ch_estim_noise_type = 'zero';
    mac_meta_rx.use_ch_estim_type = 'wiener';
    
    % channel equalization
    mac_meta_rx.use_equalization = true;

    % create rx
    rxx = dect_rx(verbose, mac_meta_rx);

    % how many antennas do we have?
    N_TX = txx.phy_4_5.tm_mode.N_TX;
    N_RX = rxx.mac_meta.N_RX;

    % let's say tx wants to send a packet: How many bits have to be provided?
    N_TB_bits = txx.phy_4_5.N_TB_bits;
    
    % bits per symbol
    N_bps = txx.phy_4_5.mcs.N_bps;
    
    % the snr might be different for each MCS
    snr_db_vec = snr_db_vec_global(cnt,:);

    for i=1:numel(snr_db_vec)
        
        current_snr = snr_db_vec(i);
        
        % extract variables for wiener weights
        if strcmp(rxx.mac_meta.use_ch_estim_type,'wiener') == true
            physical_resource_mapping_DRS_cell  = txx.phy_4_5.physical_resource_mapping_DRS_cell;
            N_b_DFT                             = txx.phy_4_5.numerology.N_b_DFT;
            N_PACKET_symb                       = txx.phy_4_5.N_PACKET_symb;
            N_b_CP                              = txx.phy_4_5.numerology.N_b_CP;
            samp_rate                           = txx.phy_4_5.numerology.B_u_b_DFT;
            noise_estim_wiener               	= 1/10^(current_snr/10);
            f_d_hertz_wiener                    = 20;
            tau_rms_sec_wiener                  = 363e-9;

            % the wiener filter for each snr
            rxx.wiener = lib_rx.channel_estimation_wiener_weights(physical_resource_mapping_DRS_cell,...
                                                                    N_b_DFT,...
                                                                    N_PACKET_symb,...
                                                                    N_b_CP, samp_rate,...
                                                                    noise_estim_wiener, f_d_hertz_wiener, tau_rms_sec_wiener);
        end

        STO_hist_coarse_local_convergence = zeros(n_bins, 1);
        STO_hist_fine_local_convergence = zeros(n_bins, 1);

        n_bits_PDC_sent_local_convergence = 0;
        n_bits_PDC_error_local_convergence = 0;
        n_packets_PDC_sent_local_convergence = 0;
        n_packets_PDC_error_local_convergence = 0;

        n_bits_PCC_sent_local_convergence = 0;
        n_bits_PCC_error_local_convergence = 0;
        n_packets_PCC_sent_local_convergence = 0;
        n_packets_PCC_error_local_convergence = 0;
    
        PER_has_converged = false;
        
        lastest_values = zeros(1,n_compare);
                
        while PER_has_converged == false
    
            % Local variables for a single MCS, required for parfor.
            % Each worker can write into these arrays.

            STO_hist_coarse_local = zeros(n_bins, n_worker);
            STO_hist_fine_local = zeros(n_bins, n_worker);

            n_bits_PCC_sent_local = zeros(1, n_worker);         % PCC
            n_bits_PCC_error_local = zeros(1, n_worker);
            n_packets_PCC_sent_local = zeros(1, n_worker);
            n_packets_PCC_error_local = zeros(1, n_worker);
            
            n_bits_PDC_sent_local = zeros(1, n_worker);         % PDC
            n_bits_PDC_error_local = zeros(1, n_worker);
            n_packets_PDC_sent_local = zeros(1, n_worker);
            n_packets_PDC_error_local = zeros(1, n_worker);

            %for j_w=1:1:n_worker
            parfor j_w=1:n_worker

                warning('off');

                % copy handle objects, changes within parfor are not permanent!
                tx = txx;
                rx = rxx;

                STO_hist_coarse = zeros(n_bins, 1);
                STO_hist_fine = zeros(n_bins, 1);

                n_bits_PCC_sent = 0;
                n_bits_PCC_error = 0;
                n_packets_PCC_error = 0;

                n_bits_PDC_sent = 0;
                n_bits_PDC_error = 0;
                n_packets_PDC_error = 0;

                % create channel
                ch                      = lib_rf_channel.rf_channel();
                ch.verbose              = verbose;
                ch.verbose_cp           = tx.phy_4_5.numerology.N_b_CP*tx.mac_meta.oversampling;
                ch.type                 = 'rayleigh';
                ch.amp                  = 1.0;
                ch.noise                = true;
                ch.snr_db             	= current_snr;
                ch.spectrum_occupied    = tx.phy_4_5.n_spectrum_occupied/tx.mac_meta.oversampling;
                ch.N_TX                	= N_TX;
                ch.N_RX               	= N_RX;
                ch.awgn_random_source   = 'global';
                ch.awgn_randomstream 	= RandStream('mt19937ar','Seed', randi(1e9,[1 1]));
                
                % random channel values (overwritten below)
                ch.d_sto                = 0;
                ch.d_cfo                = 0;
                ch.d_err_phase          = 0;
                
                ch.r_random_source      = 'global';
                ch.r_seed    	        = randi(1e9,[1 1]);

                % random channel values (overwritten below)
                ch.r_sto                = 0;
                ch.r_cfo                = 0;
                ch.r_err_phase          = 0;

                ch.r_samp_rate        	= tx.phy_4_5.numerology.B_u_b_DFT*tx.mac_meta.oversampling;
                ch.r_max_doppler     	= 19.458;                           % 1.946 19.458
                ch.r_type   	        = 'TDL-iii';
                ch.r_DS_desired         = 10^(-7.03 + 0.00*randn(1,1));
                ch.r_K                  = db2pow(9.0 + 0.00*randn(1,1));    % 93e-9;
                ch.r_interpolation      = true;
                ch.r_gains_active 	    = true;
                ch.init_rayleigh_rician_channel();

                for j=1:1:n_packets_per_worker_call
                    
                    rx.tx_handle = tx;
                    rx.ch_handle = ch;                    

                    % generate random PCC bits
                    PCC_user_bits = [];
                    if tx.mac_meta.PLCF_type == 1
                        PCC_user_bits = randi([0 1], 40, 1);
                    elseif tx.mac_meta.PLCF_type == 2
                        PCC_user_bits = randi([0 1], 80, 1);
                    end

                    % generate bits
                    PDC_user_bits = randi([0 1], N_TB_bits, 1);

                    % otherwise matlab complains about these variables not being intialized
                    PCC_user_bits_recovered = [];
                    PDC_user_bits_recovered = [];

                    % harq abort conditions
                    pcc_decoded_successfully = false;
                    pdc_decoded_successfully = false;

                    for z=0:1:max_harq_retransmissions

                        % there is a specific order for the redundany version
                        if mod(z,4) == 0
                            tx.mac_meta.rv = 0; % initial transmission
                            rx.mac_meta.rv = 0; % initial transmission
                        elseif mod(z,4) == 1
                            tx.mac_meta.rv = 2;
                            rx.mac_meta.rv = 2;
                        elseif mod(z,4) == 2
                            tx.mac_meta.rv = 3;
                            rx.mac_meta.rv = 3;
                        elseif mod(z,4) == 3
                            tx.mac_meta.rv = 1;
                            rx.mac_meta.rv = 1;
                        end

                        % let tx create the packet
                        samples_antenna_tx = tx.generate_packet(PCC_user_bits, PDC_user_bits);

                        % random STO after resampling at tx and before resamples at rx
                        sto_actual = ceil( (10+2*rand(1,1)) * rx.mac_meta.sto_config.n_samples_STF_b_os * L/M);

                        % random fractional + integer CFO
                        cfo_actual      = lib_util.pm_rand(rx.mac_meta.cfo_config.CFO_max_deviation_subcarrier_spacings);
                        cfo_actual_norm = cfo_actual * (1/(tx.phy_4_5.numerology.N_b_DFT*tx.mac_meta.oversampling)) * M/L;

                        % random phase
                        phase_actual = deg2rad( lib_util.pm_rand(360) );
                        
                        % re-randomize channel
                        ch.d_sto        = sto_actual;
                        ch.d_cfo        = cfo_actual_norm;
                        ch.d_err_phase  = phase_actual;

                        ch.r_sto        = sto_actual;
                        ch.r_cfo        = cfo_actual_norm;
                        ch.r_err_phase  = phase_actual;

                        % resample to USRP samp rate
                        if use_resampling == false

                            % we don't resample, directly pass samples through channel
                            samples_antenna_rx = ch.pass_samples(samples_antenna_tx, 0);

                            if use_ADC_scaling == true

                                % scale into ADC range as good as possible
                                samples_antenna_rx = samples_antenna_rx * ADC_range_scalingFactor;

                                % ADC
                                samples_antenna_rx = lib_rx.adc_clipping(samples_antenna_rx, n_ADC_bits);
                            end

                        else

                            samples_antenna_tx_USRP = lib_rx.resampling_polyphase(samples_antenna_tx, L, M, fir_coef_tx);

                            % pass samples through channel
                            samples_antenna_rx_USRP = ch.pass_samples(samples_antenna_tx_USRP, 0);

                            % the received antenna signal is quantizied
                            if use_ADC_scaling == true

                                % scale into ADC range as good as possible
                                samples_antenna_rx_USRP = samples_antenna_rx_USRP * ADC_range_scalingFactor;

                                % ADC
                                samples_antenna_rx_USRP = lib_rx.adc_clipping(samples_antenna_rx_USRP, n_ADC_bits);
                            end

                            % resample to DECT-2020 NR sample rate, L and M switch, also incorporates low pass filter against thermal noise at low SNRs
                            samples_antenna_rx = lib_rx.resampling_polyphase(samples_antenna_rx_USRP, M, L, fir_coef_rx);
                        end

                        % make next channel impulse response independent from this one
                        ch.reset_random_rayleigh_rician();

                        % Now let rx decode the frame.
                        % Rx can do so because it's mac_meta is the exact same.
                        [PCC_user_bits_recovered, PDC_user_bits_recovered] = rx.demod_decode_packet(samples_antenna_rx);

                        % measure the STO sync performance, but of course only if the synchronization is turned on
                        if isfield(rx.packet_data,'STO_report')
                            sto_DECT_domain = floor(sto_actual*M/L) + 1;
                            sto_coarse_error = min(rx.packet_data.STO_report.max_idx_coarse - sto_DECT_domain);
                            sto_fine_error = min(rx.packet_data.STO_report.max_idx_fine - sto_DECT_domain);
                            STO_hist_coarse = STO_hist_coarse + (histcounts(sto_coarse_error, edges))';
                            STO_hist_fine = STO_hist_fine + (histcounts(sto_fine_error, edges))';
                        end

                        % measure the BER uncoded

                        n_bits_PCC_sent = n_bits_PCC_sent + numel(tx.packet_data.pcc_enc_dbg.d);
                        n_bits_PCC_error = n_bits_PCC_error + sum(abs(double(tx.packet_data.pcc_enc_dbg.d) - double(rx.packet_data.pcc_dec_dbg.d_hard)));                

                        n_bits_PDC_sent = n_bits_PDC_sent + numel(tx.packet_data.pdc_enc_dbg.d);
                        n_bits_PDC_error = n_bits_PDC_error + sum(abs(double(tx.packet_data.pdc_enc_dbg.d) - double(rx.packet_data.pdc_dec_dbg.d_hard)));

                        % we might be done
                        if numel(PCC_user_bits_recovered) ~= 0
                            pcc_decoded_successfully = true;
                        end                

                        % we might be done
                        if numel(PDC_user_bits_recovered) ~= 0
                            pdc_decoded_successfully = true;
                        else
                            1;
                        end

                        % we continue sending retransmissions as long as not both were decoded correctly
                        if pcc_decoded_successfully == true && pdc_decoded_successfully == true
                            break;
                        end
                    end

                    % delete harq buffer
                    rx.harq_buf_40 = [];
                    rx.harq_buf_80 = [];
                    rx.harq_buf = [];

                    % check if PCC was decoded correctly, maybe there's still an error despite all the harq iterations
                    if pcc_decoded_successfully == false
                        n_packets_PCC_error = n_packets_PCC_error + 1;
                    end            

                    % check if PDC was decoded correctly, maybe there's still an error despite all the harq iterations
                    if pdc_decoded_successfully == false
                        n_packets_PDC_error = n_packets_PDC_error + 1;
                    end
                end

                % delete channel object
                delete(ch);

                STO_hist_coarse_local(:,j_w) = STO_hist_coarse;
                STO_hist_fine_local(:,j_w) = STO_hist_fine;

                n_bits_PCC_sent_local(1,j_w) = n_bits_PCC_sent;
                n_bits_PCC_error_local(1,j_w) = n_bits_PCC_error;
                n_packets_PCC_sent_local(1,j_w) = n_packets_per_worker_call;
                n_packets_PCC_error_local(1,j_w) = n_packets_PCC_error;         

                n_bits_PDC_sent_local(1,j_w) = n_bits_PDC_sent;
                n_bits_PDC_error_local(1,j_w) = n_bits_PDC_error;
                n_packets_PDC_sent_local(1,j_w) = n_packets_per_worker_call;
                n_packets_PDC_error_local(1,j_w) = n_packets_PDC_error;
            end

            STO_hist_coarse_local_convergence = STO_hist_coarse_local_convergence + sum(STO_hist_coarse_local,2);
            STO_hist_fine_local_convergence = STO_hist_fine_local_convergence + sum(STO_hist_fine_local,2);
            
            % get the latest PDC that was calculated
            n_bits_PCC_sent_local_convergence = n_bits_PCC_sent_local_convergence + sum(n_bits_PCC_sent_local);
            n_bits_PCC_error_local_convergence = n_bits_PCC_error_local_convergence + sum(n_bits_PCC_error_local);
            n_packets_PCC_sent_local_convergence = n_packets_PCC_sent_local_convergence + sum(n_packets_PCC_sent_local);
            n_packets_PCC_error_local_convergence = n_packets_PCC_error_local_convergence + sum(n_packets_PCC_error_local);
            
            % get the latest PDC that was calculated
            n_bits_PDC_sent_local_convergence = n_bits_PDC_sent_local_convergence + sum(n_bits_PDC_sent_local);
            n_bits_PDC_error_local_convergence = n_bits_PDC_error_local_convergence + sum(n_bits_PDC_error_local);
            n_packets_PDC_sent_local_convergence = n_packets_PDC_sent_local_convergence + sum(n_packets_PDC_sent_local);
            n_packets_PDC_error_local_convergence = n_packets_PDC_error_local_convergence + sum(n_packets_PDC_error_local);
                        
            % we limit the number of pakets per snr
            if n_packets_PDC_sent_local_convergence > n_max_packets_snr
                
                % assume convergence
                PER_has_converged = true;
                
            % if packet limit not reached, we still might have achieved convergence
            else
                
                % check if convergence is achieved

                % first get latest PER
                PER_latest = n_packets_PDC_error_local_convergence/n_packets_PDC_sent_local_convergence;

                % shift and add new values
                lastest_values(1:end-1) = lastest_values(2:end);
                lastest_values(end) = PER_latest;
                
                % all PERs must be above zero and for each PER a minimum number of frames must be corrupted
                if min(lastest_values) > 0 && n_packets_PDC_error_local_convergence >= n_min_corruped_packets

                    % hope for he best...
                    PER_has_converged = true;

                    % ... but check
                    for ii = 1:1:(n_compare-1)

                        a = max(lastest_values(ii), lastest_values(ii+1));
                        b = min(lastest_values(ii), lastest_values(ii+1));

                        relative_difference = abs(1-b/a);

                        if relative_difference > compare_convergence_threshold
                            PER_has_converged = false;
                        end
                    end
                end
            end
        end

        STO_hist_coarse_global(:,i) = STO_hist_coarse_local_convergence;
        STO_hist_fine_global(:,i) = STO_hist_fine_local_convergence;
        
        n_bits_PCC_sent_global(cnt,i) = n_bits_PCC_sent_local_convergence;
        n_bits_PCC_error_global(cnt,i) = n_bits_PCC_error_local_convergence;
        n_packets_PCC_sent_global(cnt,i) = n_packets_PCC_sent_local_convergence;
        n_packets_PCC_error_global(cnt,i) = n_packets_PCC_error_local_convergence;

        n_bits_PDC_sent_global(cnt,i) = n_bits_PDC_sent_local_convergence;
        n_bits_PDC_error_global(cnt,i) = n_bits_PDC_error_local_convergence;
        n_packets_PDC_sent_global(cnt,i) = n_packets_PDC_sent_local_convergence;
        n_packets_PDC_error_global(cnt,i) = n_packets_PDC_error_local_convergence;
        
        bps_global(cnt) = N_bps;
        tbs_global(cnt) = N_TB_bits;
        
        % check what the latest PER was?
        PER_latest = n_packets_PDC_error_local_convergence/n_packets_PDC_sent_local_convergence;
        
        fprintf('Convergence! MCS %d, SNR = %d , PER = %d, Packets = %d, Error Packets = %d at %s\n',...
                    cnt,...
                    current_snr,...
                    PER_latest,...
                    n_packets_PDC_sent_local_convergence,...
                    n_packets_PDC_error_local_convergence,...
                    datestr(now,'HH:MM:SS'));
        
        % if below threshold, abort entire simulation
        if PER_latest < per_abort_threshold
            break;
        end
    end
    
    fprintf('Done! MCS %d of %d at %s\n', cnt, numel(mcs_index_vec), datestr(now,'HH:MM:SS'));

    cnt = cnt + 1;
end

% check if we have old results
if isfile('results/var_only_counter.mat')
    
    % load copies
    load('results/var_only_counter.mat');
    
    % add results of old run

    STO_hist_coarse_global = STO_hist_coarse_global + STO_hist_coarse_global_cpy;
    STO_hist_fine_global = STO_hist_fine_global + STO_hist_fine_global_cpy;
    
    n_bits_PCC_sent_global = n_bits_PCC_sent_global + n_bits_PCC_sent_global_cpy;
    n_bits_PCC_error_global = n_bits_PCC_error_global + n_bits_PCC_error_global_cpy;
    n_packets_PCC_sent_global = n_packets_PCC_sent_global + n_packets_PCC_sent_global_cpy;
    n_packets_PCC_error_global = n_packets_PCC_error_global + n_packets_PCC_error_global_cpy;
    
    n_bits_PDC_sent_global = n_bits_PDC_sent_global + n_bits_PDC_sent_global_cpy;
    n_bits_PDC_error_global = n_bits_PDC_error_global + n_bits_PDC_error_global_cpy;
    n_packets_PDC_sent_global = n_packets_PDC_sent_global + n_packets_PDC_sent_global_cpy;
    n_packets_PDC_error_global = n_packets_PDC_error_global + n_packets_PDC_error_global_cpy;
end

% save results of all runs combined

STO_hist_coarse_global_cpy = STO_hist_coarse_global;
STO_hist_fine_global_cpy = STO_hist_fine_global;

n_bits_PCC_sent_global_cpy = n_bits_PCC_sent_global;
n_bits_PCC_error_global_cpy = n_bits_PCC_error_global;
n_packets_PCC_sent_global_cpy = n_packets_PCC_sent_global;
n_packets_PCC_error_global_cpy = n_packets_PCC_error_global;

n_bits_PDC_sent_global_cpy = n_bits_PDC_sent_global;
n_bits_PDC_error_global_cpy = n_bits_PDC_error_global;
n_packets_PDC_sent_global_cpy = n_packets_PDC_sent_global;
n_packets_PDC_error_global_cpy = n_packets_PDC_error_global;

save('results/var_only_counter.mat',...
        'STO_hist_coarse_global_cpy',...
        'STO_hist_fine_global_cpy',...
        'n_bits_PCC_sent_global_cpy',...
        'n_bits_PCC_error_global_cpy',...
        'n_packets_PCC_sent_global_cpy',...
        'n_packets_PCC_error_global_cpy',...
        'n_bits_PDC_sent_global_cpy',...
        'n_bits_PDC_error_global_cpy',...
        'n_packets_PDC_sent_global_cpy',...
        'n_packets_PDC_error_global_cpy');

% save all variables
save('results/var_all.mat');

%profile viewer
%profile off

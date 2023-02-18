function [json_struct, mac_meta_tx, tx, samples_antenna_tx, resampling_param] = packets_TX(filefolder, filename)

%   In Matlab, non-standard compliant features were inserted, which might have to be deactivated for a comporison against C++.
%
%   1) If symbols with STF show errors, go to file +lib_6_generic_procedures/ofdm_signal_generation_Cyclic_prefix_insertion.m and deactivate this line:
%
%       samples_antenna(1:STF_in_samples,:) = sqrt(N_eff_TX)*samples_antenna(1:STF_in_samples,:);
%
%   2) If symbols with DRS show errors, go to file +lib_5_physical_layer_transmissions/DRS.m and deactivate the following line:
%
%       y_b1 = y_b1*sqrt(N_eff_TX);
%
%   3) If symbols with PCC show errors, go to file +lib_6_generic_procedures/Transmit_diversity_precoding_Y.m and actiave this line
%
%       dont_use_preafactor_as_in_TS = true;
%
%      and deactivate this line:
%
%       dont_use_preafactor_as_in_TS = false;
%

    ffn = fullfile(filefolder, filename);
    json_struct = lib_review.load_json(ffn);

    % these are all the information we need to generate a TX packet in Matlab and compare the IQ samples with c++ result
    mac_meta_tx.u                   = json_struct.dect.u;
    mac_meta_tx.b                   = json_struct.dect.b;
    mac_meta_tx.PacketLengthType    = json_struct.dect.PacketLengthType;
    mac_meta_tx.PacketLength        = json_struct.dect.PacketLength;
    mac_meta_tx.tm_mode_0_to_11     = json_struct.dect.tm_mode;
    mac_meta_tx.mcs_index           = json_struct.dect.mcs;
    mac_meta_tx.Z                   = json_struct.dect.Z;
    mac_meta_tx.oversampling        = json_struct.dect.oversampling;
    mac_meta_tx.codebook_index      = json_struct.dect.codebook_index;
    mac_meta_tx.PLCF_type           = json_struct.dect.PLCF_type;
    mac_meta_tx.rv                  = json_struct.dect.rv;
    mac_meta_tx.network_id          = de2bi(json_struct.dect.network_id ,32,'left-msb');

    % matlab's function lteRateMatchTurbo() fails for BPSK modulation with large transport blocks, unknown why
    if mac_meta_tx.mcs_index == 0
        fprintf("\tMSC index is %f.\n", mac_meta_tx.mcs_index);
    end

    % read PLCF and TB bits from c++ (user bits pre channel coding)
    PLCF_bits_cpp = json_struct.dect.data.binary.PLCF;
    TB_bits_cpp = json_struct.dect.data.binary.TB;

    % read PCC and PDC bits from c++ (bits post channel coding)
    PCC_cpp = json_struct.dect.data.binary.PCC;
    PDC_cpp = json_struct.dect.data.binary.PDC;

    % create an equivalent matlab tx
    verbose = 0;
    tx = dect_tx(verbose, mac_meta_tx);

    % try to regenerate tx samples with the user bits from c++
    try
        samples_antenna_tx = tx.generate_packet(PLCF_bits_cpp, TB_bits_cpp);
    catch
        fprintf("\tPacket generation failed, likely lteRateMatchTurbo() threw. Check error message. Returning.\n");
        return;
    end

    % save a copy before resampling for comparison
    samples_antenna_tx_orig = samples_antenna_tx;

    % compare PCC bit errors after coding
    err_pcc_vec = abs(PCC_cpp - double(tx.packet_data.pcc_enc_dbg.d));
    err_pcc = sum(err_pcc_vec);
    if(err_pcc > 0)
        error("File with filename %s failed at PCC comparison.", filename);
    end

    % compare PDC bit errors after coding
    err_pdc_vec = abs(PDC_cpp - double(tx.packet_data.pdc_enc_dbg.d));
    err_pdc = sum(err_pdc_vec);
    if(err_pdc > 0)
        error("File with filename %s failed at PDC comparison.", filepath_and_name);
    end

    % read tx samples from json file from c++
    samples_antenna_tx_cpp = json_struct.dect.data.IQ.real + 1i*json_struct.dect.data.IQ.imag;

    % json concatenates antenna samples, reshape into matrix of size N_TX x N_samples_per_antenna
    N_samples_per_antenna = numel(samples_antenna_tx_cpp)/tx.phy_4_5.tm_mode.N_TX;
    samples_antenna_tx_cpp = reshape(samples_antenna_tx_cpp, N_samples_per_antenna, tx.phy_4_5.tm_mode.N_TX);

    % resampling?
    if json_struct.dect.resampling.L == 1 && json_struct.dect.resampling.M == 1

        resampling_param.use_resampling = false;

        % hw sample rate and DECT sample rate are the same
        hw_samp_rate                = json_struct.dect.resampling.samp_rate;
        dect_samp_rate              = hw_samp_rate;

    else

        resampling_param.use_resampling = true;

        L                           = json_struct.dect.resampling.L;
        M                           = json_struct.dect.resampling.M;
        hw_samp_rate                = json_struct.dect.resampling.samp_rate;
        normalized_input_passband   = json_struct.dect.resampling.normalized_input_passband;
        oversampling_minimum        = json_struct.dect.resampling.oversampling_minimum;

        % what is the DECT2020-NR sample rate before resampling?
        dect_samp_rate              = hw_samp_rate/L*M;

        % design resampling filter
        u = 1;
        b = 1;
        N_b_DFT = 64;

        % determine number of required guards to get the same normalized_input_passband
        n_guards_top = 32 - 2*32*normalized_input_passband*oversampling_minimum;

        [fir_coef_tx, fir_coef_rx] = lib_rx.resampling_filter_design(u, b, N_b_DFT, n_guards_top, oversampling_minimum, L, M);

        % set necessary output
        resampling_param.L = L;
        resampling_param.M = M;
        resampling_param.fir_coef_rx = fir_coef_rx;

        % resample
        samples_antenna_tx = lib_rx.resampling_polyphase(samples_antenna_tx, L, M, fir_coef_tx);

        % sanity check: make sure both files are of the same length
        len_cpp = numel(samples_antenna_tx_cpp(:,1));
        len_mat = numel(samples_antenna_tx(:,1));
        if len_cpp < len_mat
            error("C++ packet shorter than matlab packet.");
        end
        len_diff = len_cpp - len_mat;
        if abs(len_cpp - len_mat) > L
            error("Packet length difference to large.");
        end

        % adjust length of matlab packet
        samples_antenna_tx = [samples_antenna_tx; zeros(len_diff, tx.phy_4_5.tm_mode.N_TX)];
    end

    % for plotting, create two time axis (can be the same is not resampling is used)
    t_dect = (0:1:  (numel(samples_antenna_tx_orig(:,1)) - 1) ) * 1/dect_samp_rate;
    t_hw   = (0:1:  (numel( samples_antenna_tx_cpp(:,1)) - 1) ) * 1/hw_samp_rate;

    % plot comparison
    screen_size_px = get(0,'screensize');
    screen_size_x = screen_size_px(3);
    screen_size_y = screen_size_px(4);
    x_margin = 50;
    y_margin = 50;
    for i=1:1:tx.phy_4_5.tm_mode.N_TX

        fig = figure(i);

        % figure position
        if i <= 4
            fig_pos_x = x_margin + (screen_size_x - 2*x_margin) / 4 * (i-1);
            fig_pos_y = -y_margin;
            fig_pos = [fig_pos_x fig_pos_y];
        else
            fig_pos_x = x_margin + (screen_size_x - 2*x_margin) / 4 * (i-5);
            fig_pos_y = -(y_margin + screen_size_y / 2);
            fig_pos = [fig_pos_x fig_pos_y];
        end

        movegui(fig,fig_pos);
        clf()

        % compare real part
        subplot(4,1,1)
        plot(real(samples_antenna_tx_cpp(:,i)));
        hold on
        plot(real(samples_antenna_tx(:,i)), 'r');
        legend('c++','matlab');

        t = title(strcat("File Name: ", filename));
        set(t,'Interpreter','none');

        % compare imag part
        subplot(4,1,2)
        plot(imag(samples_antenna_tx_cpp(:,i)));
        hold on
        plot(imag(samples_antenna_tx(:,i)), 'r');
        legend('c++','matlab');

        % compare real before and after resampling
        subplot(4,1,3)
        plot(t_hw,   real(samples_antenna_tx_cpp(:,i)));
        hold on
        plot(t_dect, real(samples_antenna_tx_orig(:,i)), 'r');
        legend('c++ resampled','matlab before resampling');

        % compare imag before and after resampling
        subplot(4,1,4)
        plot(t_hw,   imag(samples_antenna_tx_cpp(:,i)));
        hold on
        plot(t_dect, imag(samples_antenna_tx_orig(:,i)), 'r');
        legend('c++ resampled','matlab before resampling');
    end

    % numerical samples comparison
    diff = samples_antenna_tx_cpp - samples_antenna_tx;
    diff_abs = abs(diff);
    err = max(diff_abs, [], 'all');
    if(err > 0.001)

        % At this point the samples are incorrect, so there is a bug.
        % To debug, we can take the IQ samples from C++ and try to decode the frame.
        % It will likely fail, but we can still find unexpected number in the time-frequency grid.
        % Depending on the exact TX configuration, the internals of this function must be adjusted (e.g. number of RX antennas).

        error("File with filename %s failed at samples comparison.", filepath_and_name);
    end
end

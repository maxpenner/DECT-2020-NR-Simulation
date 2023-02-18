function [antenna_streams_mapped_rev, CFO_report_integer] = sync_CFO_integer(   antenna_streams_mapped_rev_with_os_carriers,...
                                                                                STO_templates_freq_domain,...
                                                                                N_RX,...
                                                                                N_eff_TX,...
                                                                                N_b_DFT,...
                                                                                oversampling,...
                                                                                cfo_config)

    %% extract the correct STF

    switch N_eff_TX
        case 1
            STF_values_freq_domain = STO_templates_freq_domain{1};
        case 2
            STF_values_freq_domain = STO_templates_freq_domain{2};
        case 4
            STF_values_freq_domain = STO_templates_freq_domain{3};
        case 8
            STF_values_freq_domain = STO_templates_freq_domain{4};
        otherwise
            error('Unknow N_eff_TX.');
    end

    %% we need to expand the STF template if oversampling is used

    % zeros at spectrum edge
    n_oversampling_zeros = (oversampling-1)*N_b_DFT/2;

    % when oversampling, our signal must remain in the center of the spectrum
    os_idx_start = n_oversampling_zeros + 1;
    os_idx_end = os_idx_start + N_b_DFT - 1;

    STF_values_freq_domain = [zeros(n_oversampling_zeros,1); STF_values_freq_domain; zeros(n_oversampling_zeros,1)];

    %% now we have to test for each antenna, which integer CFO makes most sense

    % one metric value per possible integer CFO per RX antenna
    metric = zeros(numel(cfo_config.integer.candidate_values), N_RX);

    % go over each rx antenna
    for i=1:1:N_RX

        % extract received frequency domain samples for RX antenna i
        transmit_streams_rev_i = cell2mat(antenna_streams_mapped_rev_with_os_carriers(i));

        % we only need the STF symbol
        STF_recv = transmit_streams_rev_i(:,1);

        % Try for each possible frequency offset.
        % In Matlab, when we circshift with a positiv value, we push subcarriers "down" or towards higher indices.
        % For DECT-2020, this corrsponds to a negative integer CFO.
        % Our matrices use the same layout as DECT-2020, which is why we also experience a negative integer CFO.
        idx = 1;
        for cfo_candidate = cfo_config.integer.candidate_values

            % shift preamble
            % when cfo_candidate is negative, we assume a negative integer CFO, so we have to push our STF "down" (invert cfo_candidate)
            % when cfo_candidate is positive, we assume a positive integer CFO, so we have to push our STF "up" (invert cfo_candidate)
            STF_values_freq_domain_shifted = circshift(STF_values_freq_domain, -cfo_candidate);

            % calculate metrix
            metric(idx, i) = sum(abs(STF_values_freq_domain_shifted .* conj(STF_recv)));

            idx = idx+1;
        end
    end

    %% 

    %% combine the results from all antennas and determine the most likely integer CFO

    metric_abs = abs(metric);
    metric = sum(metric_abs,2);

    % find index of best fitting candidate
    [~, CFO_report_integer_index] = max(metric);

    % extract corresponding CFO value
    CFO_report_integer = cfo_config.integer.candidate_values(CFO_report_integer_index);

    %% the integer CFO causes a phase rotation between the OFDM symbols, which we need to determine

    % number of symbols in the packet
    N_PACKET_symb = size(transmit_streams_rev_i,2);

    % The equation for the symbol to symbol phase rotation is: 2*pi*(1 + 8/64)*integer_CFO_in_muliples_of_subcarriers.
    % 8/64 is the ratio between cyclic prefic and symbol lenght.
    sym2sym_rot = 2*pi*(1+8/64)*CFO_report_integer;

    % this line can be used for debugging to see by how much the symbols are rorated
    %sym2sym_rot_deg = rad2deg(sym2sym_rot);

    % generate rotation for entire frame
    phase_derotation = 0 : 1 : (N_PACKET_symb-1);
    phase_derotation = phase_derotation*sym2sym_rot;
    phase_derotation = exp(-1i*phase_derotation);
    phase_derotation = repmat(phase_derotation, N_b_DFT, 1);

    % this line can be used for debugging to see by how much the symbols are rorated
    %phase_derotation_deg = rad2deg(angle(phase_derotation));

    %% correct the CFO, remove os carriers, derotate phase and write into new container

    % output container without os carriers
    antenna_streams_mapped_rev = cell(N_RX,1);

    % go over each rx antenna
    for i=1:1:N_RX

        % extract received frequency domain samples for RX antenna i
        transmit_streams_rev_i = cell2mat(antenna_streams_mapped_rev_with_os_carriers(i));

        % is CFO_report_integer is negative, we need to shift everthing "up" to correct the CFO
        % is CFO_report_integer is postitve, we need to shift everthing "down" to correct the CFO
        transmit_streams_rev_i = circshift(transmit_streams_rev_i, CFO_report_integer, 1);

        % remove oversampling subcarriers
        transmit_streams_rev_i = transmit_streams_rev_i(os_idx_start : os_idx_end, :);

        % this line can be used for debugging to see by how much the symbols are rorated
        %temp0 = tx_handle.packet_data.antenna_streams_mapped{1}./transmit_streams_rev_i;
        %temp1 = rad2deg(angle(temp0));

        % derotate symbols
        transmit_streams_rev_i = transmit_streams_rev_i.*phase_derotation;

        % this line can be used for debugging to see by how much the symbols are rorated
        %temp2 = tx_handle.packet_data.antenna_streams_mapped{1}./transmit_streams_rev_i;
        %temp3 = rad2deg(angle(temp2));

        % overwrite
        antenna_streams_mapped_rev(i) = {transmit_streams_rev_i};
    end
end


function [ch_estim] = channel_estimation_perfect_SIMO(antenna_streams_mapped_rev, N_RX, N_eff_TX, tx_handle)
    
	% we need the size of the packet
    [N_b_DFT, N_PACKET_symb] = size(cell2mat(antenna_streams_mapped_rev(1)));

    % output
    ch_estim = cell(N_RX,1);

    % for each rx antenna
    for i=1:1:N_RX

        % create empty container
        ch_estim_i = zeros(N_b_DFT, N_PACKET_symb, N_eff_TX);

        % received f domain samples at antenna i
        transmit_streams_rev_i = cell2mat(antenna_streams_mapped_rev(i));

        % for each tx antenna
        for j=1:1:N_eff_TX
            %% extract send values

            % transmitted
            values = cell2mat(tx_handle.packet_data.antenna_streams_mapped(j));

            % received
            y = transmit_streams_rev_i;

            %% perfect estimation, can be NaN and +/- infinity
            ch_estim_i(:,:,j) = y./values;
        end
        ch_estim(i) = {ch_estim_i};
    end

end

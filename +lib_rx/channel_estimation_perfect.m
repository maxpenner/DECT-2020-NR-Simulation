function [ch_estim] = channel_estimation_perfect(antenna_streams_mapped_rev, N_RX, N_eff_TX, ch_handle_)
    
	% we need the size of the packet
    [N_b_DFT, N_PACKET_symb] = size(cell2mat(antenna_streams_mapped_rev(1)));

    % output
    ch_estim = cell(N_RX,1);
    
    % AWGN
    if strcmp(ch_handle_.type,'awgn') == true
        
        % for each rx antenna
        for i=1:1:N_RX

            % create empty container
            ch_estim_i = zeros(N_b_DFT, N_PACKET_symb, N_eff_TX);

            % for each tx antenna
            for j=1:1:N_eff_TX
                ch_estim_i(:,:,j) = ones(N_b_DFT, N_PACKET_symb);
            end
            ch_estim(i) = {ch_estim_i};
        end        
        
    % rayleigh
    elseif strcmp(ch_handle_.type,'rayleigh') == true || strcmp(ch_handle_.type,'rician') == true
        
        if ch_handle_.r_gains_active == false
            error('Perfect channel knowledge requires path gains.');
        end
        
        if ch_handle_.r_max_doppler ~= 0
            error('Perfect channel knowledge requires zero doppler frequency.');
        end        
        
        % only flat fading
        if numel(ch_handle_.r_matlab_MIMO_obj.AveragePathGains) == 1
            
            % for each rx antenna
            for i=1:1:N_RX

                % create empty container
                ch_estim_i = zeros(N_b_DFT, N_PACKET_symb, N_eff_TX);

                % for each tx antenna
                for j=1:1:N_eff_TX
                    ch_estim_i(:,:,j) = ch_handle_.r_gains(1,1,j,i) * ones(N_b_DFT, N_PACKET_symb);
                end
                ch_estim(i) = {ch_estim_i};
            end
            
        else
            error('Perfect channel knowledge not yet defined for non-flat fading channel type.');
        end
        
    else
        error('Channel type must be awgn, rayleigh or rician.');
    end
end

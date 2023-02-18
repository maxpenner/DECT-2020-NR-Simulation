function [ch_estim] = channel_estimation_ls(antenna_streams_mapped_rev, physical_resource_mapping_STF_cell, physical_resource_mapping_DRS_cell, N_RX, N_eff_TX)
    
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
            % linear indices of drs
            linear_indices_matlab = cell2mat(physical_resource_mapping_DRS_cell(j,4));

            % transmitted drs
            values = cell2mat(physical_resource_mapping_DRS_cell(j,3));

            % received drs
            y = transmit_streams_rev_i(ind2sub([N_b_DFT N_PACKET_symb], linear_indices_matlab));

            %% we can also include stf
            if 1==0
                % first get the position of all stf subcarriers
                STF_l_matlab    = cell2mat(physical_resource_mapping_STF_cell(2)) + MATLAB_INDEX_SHIFT;
                STF_k_i_matlab  = lib_util.index_conversion_TS_matlab(N_b_DFT, cell2mat(physical_resource_mapping_STF_cell(1)));
                STF_values      = cell2mat(physical_resource_mapping_STF_cell(3));

                % LS equalization
                ls_stf = transmit_streams_rev_i(STF_k_i_matlab,STF_l_matlab)./STF_values;

                % construct points for interpolant
                x_stf = repmat(STF_l_matlab,size(STF_k_i_matlab));
                y_stf = STF_k_i_matlab;
            else
                x_stf = [];
                y_stf = [];
                ls_stf = [];
            end

            %% least squares estimation
            ls = y./repmat(values,1, size(y,2));

            % reduce noise
            %ls = movmean(ls,[1 1],1);
            %ls = movmean(ls,[1 1],2);

            %% interpolate
        
            % rows and cols of drs values
            [y,x] = ind2sub([N_b_DFT N_PACKET_symb], linear_indices_matlab);
            x = [x_stf; reshape(x,numel(x),1)];
            y = [y_stf; reshape(y,numel(y),1)];
            ls = [ls_stf; reshape(ls,numel(ls),1)];

            % create interpolant for magnitude and angle separately
            IP_magnitude = scatteredInterpolant(x,y,abs(ls), 'linear','linear');
            IP_angle = scatteredInterpolant(x,y,ls, 'linear','linear');

            % we need interpolation at each subcarier of the packet
            [xq,yq] = meshgrid(1:1:N_PACKET_symb, 1:1:N_b_DFT);

            % interpolate
            channel_estimation_magnitude = IP_magnitude(xq, yq);
            channel_estimation_angle = angle(IP_angle(xq,yq));

            % combine
            ch_estim_i(:,:,j) = channel_estimation_magnitude.*exp(1i*channel_estimation_angle);
        end
        ch_estim(i) = {ch_estim_i};
    end
end


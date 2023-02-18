function [samples_antenna_sto, STO_CFO_report] = sync_STO(  verbose,...
                                                            u,...
                                                            N_b_DFT,...
                                                            samples_antenna,...
                                                            STF_templates_time_domain,...
                                                            n_packet_samples,...
                                                            oversampling,...
                                                            sto_config,...
                                                            cfo_config)

    %% precalculation ofparameters that are known at the receiver a-priori

    % packet size
    [n_samples_ch, N_RX] = size(samples_antenna);

    % STF pattern
    L = sto_config.n_pattern;

    % what is the length of the STF with oversampling included?
    n_samples_STF_os = sto_config.n_samples_STF_b_os;

    % sanity check
    if mod(n_samples_STF_os, L) ~= 0
        error("STF length not a multiple of L.");
    end

    % what is the length of a single pattern?
    M = n_samples_STF_os/L;

    %% coarse threshold crossing

    % For each RX antenna, we go over all samples and calculate a metric until we cross a predefined threshold.
    % We need the sample index of this theshold crossings.
    coarse_threshold_crossing_idx = zeros(N_RX, 1);

    for i=1:1:N_RX
        
        % get all samples of this particular antenna
        samples_antenna_single = samples_antenna(:,i);

        % container for detection metric
        metric = zeros(n_samples_ch - 2*n_samples_STF_os, 1);
        
        % the step must be small enough to find every STF at every SNR
        for k = 1 : sto_config.threshold.step : numel(metric)

            % determine the metric at this samples index
            metric(k) = coarse_detection(samples_antenna_single(k : k + n_samples_STF_os - 1), M, L, sto_config.minimum_power_threshold);

            % the threshold must be small enough to find every STF at every SNR
            if metric(k) >= sto_config.threshold.value

                % save where we have crossed the threshold, this is an absolute index
                coarse_threshold_crossing_idx(i) = k;

                % abort the search for this antenna
                break;
            end
        end
    end

    % keep the antennas where we found a preamble
    coarse_threshold_crossing_idx = coarse_threshold_crossing_idx(coarse_threshold_crossing_idx > 0);

    % In a real receiver, we will start the rest of the synchronization algorithm once we hit a preamble at ANY antenna.
    % This corresponds in our case to finding the smallest index across all antennas.
    % If we didn't find any preamble at all, assume index 1, which will lead to a incorrect decoding if STO is sufficiently large.
    if isempty(coarse_threshold_crossing_idx) == true
        coarse_threshold_crossing_idx = 1;
    else
        coarse_threshold_crossing_idx = min(coarse_threshold_crossing_idx);
    end

    %% after detection, extract only the range of samples what will be required for the upcoming synchronization steps

    req_start = coarse_threshold_crossing_idx + 1;
    req_end   = coarse_threshold_crossing_idx + 1 + sto_config.coarse_peak.search_length + sto_config.fine.search_area + n_samples_STF_os;

    % if we overreach, set index to 1 which will very likely to an erroneous packet (if the STO is suffiently large)
    if req_end > n_samples_ch
        req_start = 1;
        req_end   = 1 + 1 + sto_config.coarse_peak.search_length + sto_config.fine.search_area + n_samples_STF_os;
    end

    req_lenth = req_end - req_start + 1;
    
    samples_antenna_required = samples_antenna(req_start:req_end, :);

    %% each following step can be executed multiple times to increase synchonization performance
    for iter=1:1:sto_config.n_iterations

        %% coarse peak search
    
        % We now know we found an STF, now we need to find the peak of the coarse detection metric.
        % We start searching for the peak right after the index threshold_cross_idx.
        % Peaks are determined for each rx antenna individually.
    
        % how strong is the detected packet?
        coarse_peak_height = zeros(N_RX, 1);
    
        % where was the packet detected, this is now an index relative to samples_antenna_required
        coarse_peak_idx = ones(N_RX, 1);
    
        for i=1:1:N_RX
            
            % get samples of this particular antenna
            samples_antenna_single = samples_antenna_required(:,i);
    
            % container for detection metric
            metric = zeros(sto_config.coarse_peak.search_length, 1);
    
            % we continue the search with a step of one sample
            for m = 1 : sto_config.coarse_peak.step : sto_config.coarse_peak.search_length
    
                % We found the packet too late and we let the rest of the matrix at zero.
                % This can't happen in a real receiver, but here we have a limited amount of samples.
                % Necertheless, this should barely ever happen.
                if m + n_samples_STF_os-1 > n_samples_ch
                    break;
                end

                % determine the metric at this samples index
                metric(m) = coarse_detection(samples_antenna_single(m:m+n_samples_STF_os-1), M, L, sto_config.minimum_power_threshold);
            end

            % warning: if sto_config.coarse_peak.step > 1, we apply the moving average across zeros!
    
            % we apply some filtering to smooth the coarse detection
            metric_mm = movmean(metric, sto_config.coarse_peak.movmean);
    
            % finally we search for the maximum
            [cph,cpi] = max(metric_mm);
    
            % save values, these values are relative to samples_antenna_required
            coarse_peak_height(i) = cph;
            coarse_peak_idx(i)    = cpi;
    
            % debugging
            if verbose > 1
                figure()
                clf()
                hold on
    
                subplot(2,1,1)
                plot(abs(metric));
                xline(cpi, 'r');
                str = append('Coarse Synchronization Metric for antenna ', num2str(i));
                title(str);
                ylim([0 1.2]);
    
                subplot(2,1,2)
                plot(abs(metric_mm));
                xline(cpi, 'r');
                str = append('Coarse Synchronization Metric after movmean for antenna ', num2str(i));
                title(str);
                ylim([0 1.2]);
                xline(cpi, 'r');
            end
        end
    
        %% coarse peak selection
        
        % we have found a peak for each rx antenna, now we have to pick the optimal peak index
    
        % we consider only peaks above a threshold
        idx_peak_high = coarse_peak_height >= sto_config.coarse_peak.threshold;
        idx_peak_low  = coarse_peak_height < sto_config.coarse_peak.threshold;
    
        % remove peaks that are too low
        coarse_peak_height = coarse_peak_height(idx_peak_high);
        coarse_peak_idx    = coarse_peak_idx(idx_peak_high);
    
        % if can actually happen, that we have no peak above the threshold
        if isempty(coarse_peak_height) == false
    
            % the sync point are the peak indices weighted by the heiht
            normalization = coarse_peak_height/sum(coarse_peak_height);
            coarse_sync_idx_weighted = sum( coarse_peak_idx .* normalization);
            coarse_sync_idx_weighted = floor(coarse_sync_idx_weighted);
            
            % weight every peak index equally
            coarse_syc_idx_weighted_equal = floor(mean(coarse_peak_idx));
            
            % pick the peak index of the strongest peak
            [~, strongest_idx] = max(coarse_peak_height);
            coarse_sync_idx_strongest = coarse_peak_idx(strongest_idx);
    
            % we take the individial values from each antenna and fill antennas withouth detection point
            if sto_config.coarse_peak.selection_individual == true
    
                coarse_sync_idx = coarse_peak_idx;
    
                % fill antennas with peaks that are to low
                if strcmp(sto_config.coarse_peak.selection, 'weighted') == true

                    coarse_sync_idx(idx_peak_low) = coarse_sync_idx_weighted;

                elseif strcmp(sto_config.coarse_peak.selection, 'weighted equal') == true

                    coarse_sync_idx(idx_peak_low) = coarse_syc_idx_weighted_equal;

                elseif strcmp(sto_config.coarse_peak.selection, 'strongest') == true

                    coarse_sync_idx(idx_peak_low) = coarse_sync_idx_strongest;

                else
                    error('Unknown detection selection type.');
                end
    
            else
    
                % each antenna gets the same value
                if strcmp(sto_config.coarse_peak.selection, 'weighted') == true

                    coarse_sync_idx = ones(N_RX, 1)*coarse_sync_idx_weighted;

                elseif strcmp(sto_config.coarse_peak.selection, 'weighted equal') == true

                    coarse_sync_idx = ones(N_RX, 1)*coarse_syc_idx_weighted_equal;

                elseif strcmp(sto_config.coarse_peak.selection, 'strongest') == true

                    coarse_sync_idx = ones(N_RX, 1)*coarse_sync_idx_strongest;
                    
                else
                    error('Unknown detection selection type.');
                end
    
            end

            % make sure every the indices are within range
            coarse_sync_idx(coarse_sync_idx <= 0) = 1;
            coarse_sync_idx(coarse_sync_idx >= sto_config.coarse_peak.search_length) = sto_config.coarse_peak.search_length;
    
        % we didn't find any packets, assume a synchronization at 1
        else
            coarse_sync_idx = ones(N_RX, 1);
        end
    
        %% fractional cfo calculation
    
        % after coarse detection, but before fine detection, we have to determine the fractional CFO
    
        if sto_config.fine.fractional_cfo == true
    
            % create a container with the coarsly synchronized STFs only, individual for each antenna
            samples_antenna_STFs_for_fractional_cfo = zeros(n_samples_STF_os, N_RX);
            for i=1:1:N_RX
                samples_antenna_STFs_for_fractional_cfo(:,i) = samples_antenna_required(coarse_sync_idx(i) : coarse_sync_idx(i) + n_samples_STF_os - 1, i);
            end
    
            % use the STFs and determine a fractional CFO
            [~, fractional_cfo_report] = lib_rx.sync_CFO_fractional(verbose, samples_antenna_STFs_for_fractional_cfo, n_samples_STF_os, u, cfo_config);
        else
            fractional_cfo_report = zeros(N_RX, 1);
        end
    
        %% fine synchonization for each possible STF and each possible fractional + integer CFO
    
        % We have a coarse sync point, we have also determined the fractional CFO.
        % However, in order to find proper peaks for fine synchronization, we also need to correct the integer CFO.
        % The search range is limited by the crystal oscilattor precision, see +lib_rx.sync_CFO_param.m.
    
        % first create a time base which we will use to correct the fractional + integer CFO
        time_base = 0:1:(2*sto_config.fine.search_area + n_samples_STF_os - 1);
        time_base = time_base';
    
        % for integer CFO correction we need to know the search range
        n_integer_cfo_candidates = numel(cfo_config.integer.candidate_values);
    
        % N_eff_TX = {1,2,4,8}, for each we will check one stf
        metric_maxima = zeros(4, N_RX, n_integer_cfo_candidates);
        metric_maxima_index = zeros(4, N_RX, n_integer_cfo_candidates);
    
        for i=1:1:N_RX
            
            % get the samples of this particular antenna
            samples_antenna_single = samples_antenna_required(:,i);
    
            % every antenna can have it's own search range
            search_area_min = max(1, coarse_sync_idx(i) - sto_config.fine.search_area);
            search_area_max = coarse_sync_idx(i) + sto_config.fine.search_area;
    
            % lenght of search field
            n_search_area = search_area_max - search_area_min + n_samples_STF_os;
            
            % try every possible STF type (can be limited by the maximum number of antennas)
            for j=1:1:numel(STF_templates_time_domain)
    
                % stf templates are already oversampled
                STF_template_candidate = cell2mat(STF_templates_time_domain(j));
    
                % normalize power of stf
                STF_template_candidate = STF_template_candidate/rms(STF_template_candidate);
    
                % try every possible integer CFO value
                for k=1:1:n_integer_cfo_candidates
    
                    % extract the search range in which we will be looking for the fine sync peaks
                    samples_antenna_single_search_range = samples_antenna_single(search_area_min : search_area_max + n_samples_STF_os - 1, :);
    
                    % fractional CFO has one value per antenna, integer CFO is assumed to be the same
                    total_cfo = fractional_cfo_report(i) + cfo_config.integer.candidate_values(k)*1/(N_b_DFT*oversampling);
    
                    % derotate samples in search range
                    samples_antenna_single_search_range = samples_antenna_single_search_range.*exp(1i*2*pi*(-total_cfo)*time_base(1:n_search_area));
    
                    % perform a cross correlation between the samples and the stf template
                    metric = abs(lib_rx.sync_xcorr(samples_antenna_single_search_range, STF_template_candidate));
        
                    % debugging
                    if verbose > 2
                        figure()
                        plot(abs(metric));
                        str = append('Synchronization Metric for antenna ', num2str(i), ' and stf ', num2str(j));
                        title(str);
                    end
        
                    % this maximum is local
                    [~, max_idx] = max(abs(metric));
    
                    % save the maximum
                    metric_maxima(j,i,k) = abs(metric(max_idx));
    
                    % this index is converted to global
                    max_idx = max_idx + search_area_min - 1;
        
                    % save global maximum
                    metric_maxima_index(j,i,k) = max_idx;
                end
            end
        end
    
        %% find best fitting integer CFO
    
        % find the largest peak across all possible STFs and rx antennas
        [~, linear_idx] = max(metric_maxima, [], 'all', 'linear');
        [~,~,integer_cfo_best_index] = ind2sub(size(metric_maxima),linear_idx);
    
        % extract the optimal value in multiples of subcarriers
        integer_cfo_report = cfo_config.integer.candidate_values(integer_cfo_best_index);
    
        % reduce matrix size by dropping all values for all the other integer CFOs
        metric_maxima = metric_maxima(:,:,integer_cfo_best_index);
        metric_maxima_index = metric_maxima_index(:,:,integer_cfo_best_index);
        
        %% find best fitting N_eff_TX
    
        if strcmp(sto_config.fine.N_eff_TX_selection, 'weighted') == true
    
            % average peaks across antennas
            metric_maxima_sum = sum(metric_maxima,2);
            
            % find the best fitting stf
            [~, j_best] = max(metric_maxima_sum);
    
        elseif strcmp(sto_config.fine.N_eff_TX_selection, 'strongest') == true
    
            % find the maximum across all antennas
            [~, idx] = max(metric_maxima, [], 'all');
    
            [j_best, ~] = ind2sub([4 N_RX],idx);
    
        else
            error('Unknown N_eff_TX_report type.');
        end
        
        % how many effective antennas do we have?
        switch j_best
            case 1
                N_eff_TX_report = 1;
            case 2
                N_eff_TX_report = 2;
            case 3
                N_eff_TX_report = 4;
            case 4
                N_eff_TX_report = 8;
        end
    
        %% fine best fitting synchronization point
    
        metric_maxima_fine = metric_maxima(j_best, :);
    
        max_idx_fine = metric_maxima_index(j_best, :);
    
        % weight a detection point by the height of the coarse metric of the respective antenna
        if sum(metric_maxima_fine) > 0
            normalization = metric_maxima_fine/sum(metric_maxima_fine);
        else
            normalization = 1/numel(metric_maxima_fine);
        end
        max_idx_fine_weighted = sum( max_idx_fine .* normalization);
        max_idx_fine_weighted = floor(max_idx_fine_weighted);
        
        % weight every detection point equally
        max_idx_fine_weighted_equal = floor(mean(max_idx_fine));
        
        % pick the detection point with the heightest value of its coarse metric
        [~, strongest_idx] = max(metric_maxima_fine);
        max_idx_fine_strongest = max_idx_fine(strongest_idx);
    
        % calculate some metric points
    
        % fill antennas without detection
        if strcmp(sto_config.fine.selection, 'individual') == true

            % do nothing

        elseif strcmp(sto_config.fine.selection, 'weighted') == true

            max_idx_fine = ones(N_RX, 1)*max_idx_fine_weighted;

        elseif strcmp(sto_config.fine.selection, 'weighted equal') == true

            max_idx_fine = ones(N_RX, 1)*max_idx_fine_weighted_equal;

        elseif strcmp(sto_config.fine.selection, 'strongest') == true

            max_idx_fine = ones(N_RX, 1)*max_idx_fine_strongest;

        else
            error('Unknown detection selection type.');
        end

        %% if there are still interation to come, correct the CFO as we know it

        if iter < sto_config.n_iterations

            % we need a new time base with the length of one frame
            time_base = 0:1:(req_lenth - 1);
            time_base = time_base';
    
            for i=1:1:N_RX
        
                % we combine fractional and integer cfo
                total_cfo = fractional_cfo_report(i) + integer_cfo_report*1/(N_b_DFT*oversampling);
    
                % extract samples of packet and derorate
                samples_antenna_required(:,i) = samples_antenna_required(:, i) .* exp(1i*2*pi*(-total_cfo)*time_base);
            end
        end
    end

    %% convert local indices into global indices

    max_idx_coarse = coarse_sync_idx + coarse_threshold_crossing_idx;
    max_idx_fine = max_idx_fine + coarse_threshold_crossing_idx;

    %% extract one packet starting at the fine synchronization point and correct the CFO
    
    % packet is longer when oversampling
    n_packet_samples = n_packet_samples*oversampling;

    % we need a new time base with the length of one frame
    time_base = 0:1:(n_packet_samples - 1);
    time_base = time_base';
    
    % output container
    samples_antenna_sto = zeros(n_packet_samples, N_RX);
    
    % extract synchronized frames for each antenna
    for i=1:1:N_RX
        
        % minimum is 1
        if max_idx_fine(i) <= 0
            max_idx_fine(i) = 1;
        end
		
        max_idx_fine_plus_packet = max_idx_fine(i) + n_packet_samples - 1;

        % make sure we don't overreach
        if max_idx_fine_plus_packet > n_samples_ch
            max_idx_fine(i) = n_samples_ch - n_packet_samples + 1;
            max_idx_fine_plus_packet = n_samples_ch;
        end

        % extract synchronized packet
        if sto_config.cfo_correction.fractional == true || sto_config.cfo_correction.integer == true

            total_cfo = 0;

            if sto_config.cfo_correction.fractional == true
                total_cfo = fractional_cfo_report(i);
            end

            if sto_config.cfo_correction.integer == true
                total_cfo = total_cfo + integer_cfo_report*1/(N_b_DFT*oversampling);
            end

            % extract samples of packet and derorate
            samples_antenna_sto(:,i) = samples_antenna(max_idx_fine(i) : max_idx_fine_plus_packet, i) .* exp(1i*2*pi*(-total_cfo)*time_base);
        else

            % extract samples of packet
            samples_antenna_sto(:,i) = samples_antenna(max_idx_fine(i) : max_idx_fine_plus_packet, i);
        end
    end

    %% create output report
    STO_CFO_report.max_idx_coarse   = max_idx_coarse;
    STO_CFO_report.N_eff_TX         = N_eff_TX_report;
    STO_CFO_report.fractional_cfo   = fractional_cfo_report;
    STO_CFO_report.integer_cfo      = integer_cfo_report;
    STO_CFO_report.max_idx_fine     = max_idx_fine;
end

function metric = coarse_detection(window, M, L, detection_E_rms_threshold)
    
    % separate into repetitive pattern
    window = reshape(window, M, L);

    % power for normalization
    E = sum(sum(window.*conj(window)));

    % we have to have a minimal power to assume there is a frame
    E_mean = sqrt(E/(M*L));
    if E_mean < detection_E_rms_threshold
        metric = 0;
        return;
    end
    
    % correlation from pattern to pattern
    P = 0;
    for j=1:1:L-1
        P = P + sum(window(:,j).*conj(window(:,j+1)));
    end

    metric = (  (L/(L-1)) * abs(P)/E  )^2;

    % We assume the worst case:
    % We can't detect NaN and inf, instead we get maximum metric.
    if isnan(metric) || isinf(metric)
        metric = 1;
    end
end


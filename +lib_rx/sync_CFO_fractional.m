function [samples_antenna_sto_cfo, CFO_report] = sync_CFO_fractional(verbose, samples_antenna_sto, n_STF_template, u, cfo_config)

    % Technical Specification assumes first index is 0, matlab 1
    MATLAB_INDEX_SHIFT = 1;

    % SOURCES
    % Schmidl Cox: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=650240
    % BLUE: https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=752907
    % wiki: https://en.wikipedia.org/wiki/Carrier_frequency_offset
    
    % how often does the pattern repeat?
    switch u
        case 1
            L = 7;
        case {2,4,8}
            L = 9;
    end
    
    % number of samples in a single pattern repetition
    n_STF_pattern = n_STF_template/L;
    
    % get the input size of the input packet
    [n_samples_antenna_sto, N_RX] = size(samples_antenna_sto);
    
    % we have one CFO estimation per RX antenna, it does not depend on the number of TX antennas
    CFO_report = zeros(N_RX,1);

    % output container
    samples_antenna_sto_cfo = zeros(size(samples_antenna_sto));
    
    % time base for CFO correction
    time_base = 0:1:(n_samples_antenna_sto-1);
    time_base = time_base';

    % if we average the cfo estimation across multiples antennas, we require a container for the estimation metric
    if cfo_config.fractional.correct_all_antennas_individually == false

        if strcmp(cfo_config.type, 'SchmidlCox') == true
            P_total = 0;
        elseif strcmp(cfo_config.type, 'BLUE') == true
            R_m_total = 0;
        else
            error('Unknown CFO estimation type.');
        end
    end
    
    % process each antenna
    for i=1:1:N_RX
        
        % extract the preamble as we have received it, perfect symbol synchronization is assumed
        STF_found = samples_antenna_sto(1:n_STF_template, i);

        if strcmp(cfo_config.type, 'SchmidlCox') == true

            % SCHMIDL COX
            % Schmidl Cox assumes one OFDM symbol that consists of 2 equal halfes.
            % DECT-2020 has 1.75 (u=1) or 2.25 (u=2,4,8) OFDM symbols with 7 or 9 repetitions.
            % Our capture rate is therefore doubled in comparison to Schmidl Cox, we can capture a CFO of +/- two subcarrier spacings.
            % Schmidl Cox can capture only +/- one subcarrier spacing.

            % reorder to repeated patterns
            STF_found_mat = reshape(STF_found, n_STF_pattern, L);
            
            % equation (5)
            P = STF_found_mat(:,2:end).*conj(STF_found_mat(:,1:end-1));
    
            % summation equation (5)
            P = sum(sum(P));
            
            % estimate the cfo at this antenna
            if cfo_config.fractional.correct_all_antennas_individually == true
                CFO_report(i) = cfo_from_complex_samples_SchmidlCox(P, n_STF_pattern);
            else
                P_total = P_total + P;
            end

        elseif strcmp(cfo_config.type, 'BLUE') == true

            % BLUE
            % The following is very confusing at the beginning:
            %
            %   The paper on the BLUE estimator assumes that N is the number of samples in one OFDM symbols after CP removal.
            %   In our case, N is the length of 1.75 or 2.25 OFDM symbols (depending on u) with oversampling.
            %   However, we always have 4 repetitions in one OFDM symbol, so we can capture only +/- two subcarrier spacings.
    
            % BLUE algorithm as presented in https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=752907
            N = n_STF_template;     % in the paper this is one OFDM symbol after CP removal, for us this is one symbol + CP
            M = n_STF_pattern;
            H = ceil(L/2);          % H = L/2 is supposed to be optimal
    
            % equation (3)
            m_vec = 0:H;
            R_m = zeros(size(m_vec));
            for m=m_vec
                k = (m*M) : 1 : (N-1);
                A = STF_found(k + MATLAB_INDEX_SHIFT);
                B = STF_found(k - m*M + MATLAB_INDEX_SHIFT);
                C = A.*conj(B);
                D = sum(C);
                
                % UNCLEAR: What is the purpose of the factor 1/(N-m*M)? We're taking the angle(.) anyway.
                R_m(m+MATLAB_INDEX_SHIFT) = 1/(N-m*M)*D;
            end
    
            if cfo_config.fractional.correct_all_antennas_individually == true
                CFO_report(i) = cfo_from_complex_samples_BLUE(L, H, R_m, N);
            else
                R_m_total = R_m_total + R_m;
            end

        else
            error('Unknown CFO estimation type.');
        end

        % derotate carrier frequency offset for this particular antenna
        if cfo_config.fractional.correct_all_antennas_individually == true
            samples_antenna_sto_cfo(:,i) = samples_antenna_sto(:,i).*exp(1i*2*pi*(-CFO_report(i))*time_base);
        end
    end

    % derotate carrier frequency offset for each antenna with one common CFO estimation
    if cfo_config.fractional.correct_all_antennas_individually == false

        if strcmp(cfo_config.type, 'SchmidlCox') == true
            CFO_report = ones(N_RX,1) * cfo_from_complex_samples_SchmidlCox(P_total, n_STF_pattern);
        elseif strcmp(cfo_config.type, 'BLUE') == true
            CFO_report = ones(N_RX,1) * cfo_from_complex_samples_BLUE(L, H, R_m_total, N);
        else
            error('Unknown CFO estimation type.');
        end
        
        for i=1:1:N_RX
            samples_antenna_sto_cfo(:,i) = samples_antenna_sto(:,i).*exp(1i*2*pi*(-CFO_report(i))*time_base);
        end
    end
end

function cfo_report = cfo_from_complex_samples_SchmidlCox(P, n_STF_pattern)
    
    % by how much does the fractional STO rotate over one pattern normalized to 2*pi?
    cfo_report = angle(P)/(2*pi);

    % this line can be used to see the rotation over one OFDM symbol without CP
    %cfo_report_one_symbol = cfo_report*4;

    % The CFO we need to use for correction at a sampling rate of 1.
    % Note that N_b_DFT_os = 4*n_STF_pattern.
    cfo_report = cfo_report/n_STF_pattern;
end

function cfo_report = cfo_from_complex_samples_BLUE(L, H, R_m, N)

    % Technical Specification assumes first index is 0, matlab 1
    MATLAB_INDEX_SHIFT = 1;

    % equation (7)
    m_vec = 1:H;
    phi = zeros(size(m_vec));
    for m=m_vec
        A = angle(R_m(m + MATLAB_INDEX_SHIFT));
        B = angle(R_m(m-1 + MATLAB_INDEX_SHIFT));
        C = A-B;
        phi(m) = lib_util.wrap2_plus_pi_minus_pi(C);
    end
    
    % DEBUGGING
    %cfo_estimations_of_phi = phi/(2*pi*M);
    
    % equation (16)
    m_vec = 1:H;
    w_m = zeros(size(m_vec));
    for m=m_vec
        w_m(m) = 3*(  (L-m)*(L-m+1) - H*(L-H)  ) / (  H*(4*H^2 - 6*L*H + 3*L^2 - 1)  );
    end
    
    % check if sum of weights is close to one
    assert(1.0001 > sum(w_m) && sum(w_m) > 0.9999, 'Weights are not close enough to 1');
    
    % equation (13)
    % This equation is not stable for CFOs of 2, 6, 10 and 14 subcarriers (an negatives thereof).
    % In these cases, phi is very close to pi.
    % But due to noise, the individual values can be +pi or -pi and then the weighted sum can become values far off those integer CFOs.
    cfo_report = 1/(2*pi/L) * sum(w_m .* phi);
    
    % the paper on the BLUE estimator assumes a cfo normalized to subcarrier spacing, so we have to mulitply by N
    % see equation (1)
    cfo_report = cfo_report/N;
end


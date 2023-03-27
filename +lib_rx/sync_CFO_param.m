function [cfo_config] = sync_CFO_param(u)

    %% FRACTIONAL (pre FFT, based on STF)

    % SchmidlCox for an OFDM symbol with 4 repetitions instead of 2, capture range is +/- 2 subcarriers instead of +/- 1 subcarrier

    %% INTEGER (based on STF)
    %
    % How large is the search space for the integer CFO?
    %
    %   According to the DECT-2020 NR standard, battery powered devices are allowed to have up to 30ppm.
    %   At the receiver, we see a maximum of 2*30ppm = 60ppm deviation.
    %   At 6GHz, 60 ppm corresponds to 360kHz.
    %   The minimum subcarrier spacing is 27kHz, so the maximum CFO is 360kHz/27kHz = 13.33333 subcarrier spacings.
    %
    %   With larger u, the subcarrier spacing is larger as well and the actual possible deviation becomes smaller.
    
    % get the maximum physical deviation in muliples of the subcarrie spacing in use
    cfo_config.CFO_max_deviation_subcarrier_spacings = sync_CFO_max_estimation(u);
    
    % saw tooth of fractional CFO correction
    if cfo_config.CFO_max_deviation_subcarrier_spacings < 2
        search_lim = 0;
    elseif cfo_config.CFO_max_deviation_subcarrier_spacings < 6
        search_lim = 4;
    elseif cfo_config.CFO_max_deviation_subcarrier_spacings < 10
        search_lim = 8;
    elseif cfo_config.CFO_max_deviation_subcarrier_spacings < 14
        search_lim = 12;
    end
    
    % we add some more possible subcarrier deviations for testing as the fractional CFO correction is not optimal
    search_lim = search_lim + 4;
    
    % search exctly within this range
    %cfo_config.integer.candidate_values = -search_lim : 1 : search_lim;
    cfo_config.integer.candidate_values = -search_lim : 4 : search_lim;

    %% RESIDUAL (post FFT, based on STF and DRS)
    
    % nothing to define here
end

function [CFO_max_deviation_subcarrier_spacings] = sync_CFO_max_estimation(u)

    % sanity check
    if ismember(u, [1 2 4 8]) == false
        error('Unknown u.');
    end

    % subcarrier spacing
    subc_spacing = 27000 * u;
    
    % maximum operational frequency for DECT-2020 NR
    f_max = 6e9;

    % according to table 5.5.2-1 in part 2, battery powered devices can have up to 30 ppm in extreme conditions
    ppm_max = 30;
    
    % number of subcarriers that our signal can deviate (2 as we can have 30 ppm at the transmitter and 30 ppm at the receiver)
    CFO_max_deviation_subcarrier_spacings = f_max * 2 * ppm_max/1e6 / subc_spacing;
end


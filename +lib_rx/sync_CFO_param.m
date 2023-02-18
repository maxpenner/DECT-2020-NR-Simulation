function [cfo_config] = sync_CFO_param(u)

    %% FRACTIONAL (pre FFT, based on STF)

    % do we correct the CFO derived from the STF for each RX antenna individually (true) or do we average the CFO across all RX antennas (false)?
    cfo_config.fractional.correct_all_antennas_individually = false;

    % SchmidlCox for an OFDM symbol with 4 repetitions instead of 2, capture range is +/- 2 subcarriers instead of +/- 1 subcarrier
    cfo_config.type = 'SchmidlCox';

    % It's best not to use BLUE, as it isn't stable around CFOs of 2, 6, 10 and 14 subcarriers.
    % Also, the calculations are fairly complex and therefore not well controllable.
    %cfo_config.type = 'BLUE';

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
    cfo_config.CFO_max_deviation_subcarrier_spacings = lib_rx.sync_CFO_max_estimation(u);
    
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
    % TODO
end


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


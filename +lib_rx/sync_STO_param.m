function [sto_config] = sync_STO_param(u, b, oversampling)

    % what is the length of the STF in samples for u=1 and b=1 and how many pattern does it contain?
    switch u
        case 1
            n_samples_STF = (64+8) * 14/9;
            n_pattern = 7;
        case {2,4,8}
            n_samples_STF = (64+8) * 2;
            n_pattern = 9;
        otherwise
            error('Unknown u.');
    end

    % save
    sto_config.n_samples_STF = n_samples_STF;
    sto_config.n_pattern = n_pattern;

    % what is the actual size of the STF with oversampling?
    sto_config.n_samples_STF_b_os = n_samples_STF*b*oversampling;

    % based on the length of the STF we can estimate the length of the coarse detection metric without noise
    n_samples_coarse_metric_first_half_no_noise = n_samples_STF * b * oversampling * (n_pattern-1) / n_pattern;

    % these are experience values

    % one iteration works best
    sto_config.n_iterations                     = 1;

    % must be low enough to detect at very low SNRs, must be high enough to avoid numerical imprecision
    sto_config.minimum_power_threshold          = 0.001;

    sto_config.threshold.step                   = 10*b*oversampling;
    sto_config.threshold.value                  = 0.15;

    sto_config.coarse_peak.search_length        = 1.75*n_samples_coarse_metric_first_half_no_noise;
    sto_config.coarse_peak.step                 = 1;
    sto_config.coarse_peak.threshold            = 0.15;
    sto_config.coarse_peak.movmean              = [6*b*oversampling 2*b*oversampling];
    sto_config.coarse_peak.selection_individual = false;
    sto_config.coarse_peak.selection            = 'weighted'; % 'weighted' 'weighted equal' 'strongest'

    sto_config.fine.fractional_cfo              = true;
    sto_config.fine.search_area                 = 24*b*oversampling;
    sto_config.fine.N_eff_TX_selection          = 'weighted'; % 'weighted' 'strongest'
    sto_config.fine.selection                   = 'weighted'; % 'individual' 'weighted' 'weighted equal' 'strongest'

    % Internally, both fractional and integer CFO need to be corrected to enable proper synchronization.
    % However, afte fine synchronization, we can correct the fractional CFO with optimal synchronization and integer CFO in frequency domain.
    sto_config.cfo_correction.fractional        = false;
    sto_config.cfo_correction.integer           = false;
end


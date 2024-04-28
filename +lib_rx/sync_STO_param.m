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
    sto_config.n_samples_STF_cp_only_b_os = sto_config.n_samples_STF_b_os - 64*b*oversampling;

    % based on the length of the STF we can estimate the length of the coarse detection metric without noise
    n_samples_coarse_metric_first_half_no_noise = n_samples_STF * b * oversampling * (n_pattern-1) / n_pattern;

    % The following values define parameters for the STF time domain synchronization, i.e. the search for where a packet starts.
    % All values are experience values, which should for work for all DECT-2020 NR packet configurations.
    % In a real receiver, some values could be optimized for the given hardware, for instance, threshold values matching ADC resolution.
    % These variables are explained and used in +lib_rx.sync_STF.m.

    % ################################################
    % coarse metric detection, based on autocorrelation of incoming samples

    % power threshold: must be low enough to detect even at very low SNRs, but high enough to avoid numerical imprecision
    sto_config.minimum_power_threshold      = 0.001;

    % largest step is 16*b*oversampling, i.e. one STF pattern
    sto_config.threshold.step               = 8*b*oversampling;    

    % coarse metric is normalized between 0 and 1.0, threshold should be low enough to detect at low SNR, but not too low to avoid false alarms due to noise
    sto_config.threshold.value              = 0.15;

    % This parameter has become necessary with the newly introduced cover sequence.
    % The cover sequence has made the coarse metric very narrow, so it can happen that the metric is detected on the falling edge after the coarse peak.
    % As a countermeasure, we jump back some samples at the detection point. This way we make sure that the coarse peak search begin BEFORE the coarse peak.
    sto_config.threshold.jump_pack          = sto_config.n_samples_STF_b_os / n_pattern * 2;

    % ################################################
    % coarse peak search after detection, based on autocorrelation of incoming samples

    % starting from the detection point, the search length must be long enough to definitely contain the coarse peak
    sto_config.coarse_peak.search_length    = round(sto_config.threshold.jump_pack + 1.2*n_samples_coarse_metric_first_half_no_noise);

    % when using oversampling, this step can be made larger than 1
    sto_config.coarse_peak.step             = 1;
    sto_config.coarse_peak.threshold        = 0.15;

    % additional smoothing of the coarse metric
    sto_config.coarse_peak.movmean          = [3*b*oversampling 3*b*oversampling];

    % ################################################
    % fine peak search around the coarse peak, based on crosscorrelation with precalculated STF templates
    sto_config.fine.search_area             = 24*b*oversampling;
end


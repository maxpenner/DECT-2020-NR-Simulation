function [fir_coef_tx, fir_coef_rx] = resampling_filter_design(u, b, N_b_DFT, n_guards_top, oversampling, L, M)

    % At Nyquist, our passband is 0.5 Hz at a sample rate of 1 Hz.
    % However, due to guards bandswe don't use the full passband. It is slightly less than 0.5 Hz.
    Nyquist_passband_optimized = (N_b_DFT/2 - n_guards_top)/N_b_DFT;

    % sanity check
    if L <= M
        error('For resampling L must be larger than M.');
    end

    % Image filter and ani aliasing filter merge into one.
    % Since L > M, we use L at both TX and RX for determining the passband, as this creates a narrower filter.

    % TX resampling filter design

    % design a low pass filter for resampling at the transmitter
    if oversampling == 1
        
        if Nyquist_passband_optimized >= 0.5
            error('For critically sampled signals passband must end below 0.5 Hz.');
        end

        % critically sampled
        f_pass = Nyquist_passband_optimized * 1/L;
        f_stop = (1.0 - Nyquist_passband_optimized) * 1/L;

    else
        
        % Option A
        Nyquist_passband_used = Nyquist_passband_optimized;
        
        % Option B
        %Nyquist_passband_used = 0.5;

        % this is a far better solution for a low number of fir samples, oversampling should be at least 2
        f_pass = Nyquist_passband_used / oversampling * 1/L;
        f_stop = (1 - Nyquist_passband_used / oversampling) * 1/L;

    end

    % generate a filter with Matlab function
    fir_tx = designfilt('lowpassfir', ...
                        'PassbandFrequency',f_pass, ...
                        'StopbandFrequency',f_stop, ...
                        'PassbandRipple',5, ...
                        'StopbandAttenuation',30, ...
                        'DesignMethod','kaiserwin', ...
                        'SampleRate',1);

    % extract filter coefficients
    fir_coef_tx_reference = L*tf(fir_tx);

    % generate filter coefficients with custom function
    fir_coef_tx = lib_rx.resampling_kaiser( fir_tx.PassbandFrequency,...
                                            fir_tx.StopbandFrequency,...
                                            fir_tx.PassbandRipple,...
                                            fir_tx.StopbandAttenuation,...
                                            fir_tx.SampleRate,...
                                            true);

    fir_coef_tx = L*fir_coef_tx;

    % RX resampling filter design

    % design a low pass filter for resampling at the receiver
    if oversampling == 1

        if Nyquist_passband_optimized >= 0.5
            error('For critically sampled signals passband must end below 0.5 Hz.');
        end

        % critically sampled
        f_pass = Nyquist_passband_optimized * 1/L;
        f_stop = (1.0 - Nyquist_passband_optimized) * 1/L;

    else

        % Option A
        %Nyquist_passband_used = Nyquist_passband_optimized;
        
        % Option B
        Nyquist_passband_used = 0.5;

        % this is a far better solution for a low number of fir samples, oversampling should be at least 2
        f_pass = Nyquist_passband_used / oversampling * 1/L;

        %f_stop = 0.4*1/L_or_M_larger;
        f_stop = f_pass + f_pass/(u*b);
    end

    % generate a filter with Matlab function
    fir_rx = designfilt('lowpassfir', ...
                        'PassbandFrequency',f_pass, ...
                        'StopbandFrequency',f_stop, ...
                        'PassbandRipple',5, ...
                        'StopbandAttenuation',30, ...
                        'DesignMethod','kaiserwin', ...
                        'SampleRate',1);



    % extract filter coefficients
    fir_coef_rx_reference = M*tf(fir_rx);

    % generate filter coefficients with custom function
    fir_coef_rx = lib_rx.resampling_kaiser( fir_rx.PassbandFrequency,...
                                            fir_rx.StopbandFrequency,...
                                            fir_rx.PassbandRipple,...
                                            fir_rx.StopbandAttenuation,...
                                            fir_rx.SampleRate,...
                                            true);

    fir_coef_rx = M*fir_coef_rx;

    % debug
    if 1==0

        % TX filter

        % plot filter with matlab function
        fvtool(fir_tx)
    
        % plot filters im time domain for comparison
        figure()
        clf();

        % reference filter matlab
        N = numel(fir_coef_tx_reference);
        n = 0:1:(N-1);
        plot(n-(N-1)/2, fir_coef_tx_reference, 'LineWidth', 2.2)
        hold on

        % custom filter
        N = numel(fir_coef_tx);
        n = 0:1:(N-1);
        plot(n-(N-1)/2, fir_coef_tx, 'r')
        legend('Matlab Reference', 'Custom Function')
        title('TX Resampling Filter')

        % RX filter

        % plot filter with matlab function
        fvtool(fir_rx)

        % plot filters im time domain for comparison
        figure()
        clf();

        % reference filter matlab
        N = numel(fir_coef_rx_reference);
        n = 0:1:(N-1);
        plot(n-(N-1)/2, fir_coef_rx_reference, 'LineWidth', 2.2)
        hold on

        % custom filter
        N = numel(fir_coef_rx);
        n = 0:1:(N-1);
        plot(n-(N-1)/2, fir_coef_rx, 'r')
        legend('Matlab Reference', 'Custom Function')
        title('RX Resampling Filter')
    end
end


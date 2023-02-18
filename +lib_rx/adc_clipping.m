function [samples_antenna_rx_clipped] = adc_clipping(samples_antenna_rx, n_ADC_bits)

    re = real(samples_antenna_rx);
    im = imag(samples_antenna_rx);
    
    % clip signal
    re(re > 1)  =  1;
    re(re < -1) = -1;
    im(im > 1)  =  1;
    im(im < -1) = -1;
    samples_antenna_rx_clipped = re + 1i*im;

    % ADC configuration
    qLevels = 2^n_ADC_bits;
    signalMin = -1;             % SDR outputs floats from -1 to +1
    signalMax = 1;
    scalingFactor = (signalMax-signalMin)/qLevels;
    
    % now apply ADC quantization
    samples_antenna_rx_clipped = samples_antenna_rx_clipped/scalingFactor;
    samples_antenna_rx_clipped = round(samples_antenna_rx_clipped);
    samples_antenna_rx_clipped = samples_antenna_rx_clipped*scalingFactor;
    
    % debug
    if 1==0
        figure()
        clf()

        subplot(2,1,1)
        plot(real(samples_antenna_rx_clipped))
        ylim([-1.1 1.1]);

        subplot(2,1,2)
        plot(imag(samples_antenna_rx_clipped))
        ylim([-1.1 1.1]);
    end
end
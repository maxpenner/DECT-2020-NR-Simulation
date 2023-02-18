clear all;
close all;
clc;

% debug
% fir_coeff = blackman_fir(0.1,0.08);
% plot(fir_coeff);

% sources
%   https://tomroelandts.com/articles/how-to-create-a-simple-low-pass-filter
%   https://tomroelandts.com/articles/how-to-create-a-configurable-filter-using-a-kaiser-window

% kaiser definition
PassbandFrequency   = 0.0227;
StopbandFrequency   = 0.0773;
%StopbandFrequency   = 0.0450;
PassbandRipple      = 0.23;
StopbandAttenuation = 30;
SampleRate          = 1;

% create matlab reference coefficients
LpFilt = designfilt('lowpassfir', ...
                    'PassbandFrequency',PassbandFrequency*SampleRate, ...
                    'StopbandFrequency',StopbandFrequency*SampleRate, ...
                    'PassbandRipple',PassbandRipple, ...
                    'StopbandAttenuation',StopbandAttenuation, ...
                    'SampleRate',SampleRate, ...
                    'DesignMethod','kaiserwin');

%fvtool(LpFilt)
fir_coef_reference = LpFilt.Coefficients;

% create custom fir filter
fir_coef = lib_rx.resampling_kaiser(    PassbandFrequency*SampleRate,...
                                        StopbandFrequency*SampleRate,...
                                        PassbandRipple,...
                                        StopbandAttenuation,...
                                        SampleRate,...
                                        true);

% plot filter
figure(1)
clf();

N = numel(fir_coef_reference);
n = 0:1:(N-1);
plot(n-(N-1)/2, fir_coef_reference, 'LineWidth', 2.2)

hold on

N = numel(fir_coef);
n = 0:1:(N-1);
plot(n-(N-1)/2, fir_coef, 'r')

legend('Matlab Reference', 'Custom Function')

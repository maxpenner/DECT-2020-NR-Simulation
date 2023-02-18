clear all;
close all;

% design resampling filter
u = 1;
b = 1;
N_b_DFT = 64;
n_guards_top = 3;
oversampling = 2;
L = 10;
M = 9;
[fir_coef_tx, fir_coef_rx] = lib_rx.resampling_filter_design(u, b, N_b_DFT, n_guards_top, oversampling, L, M);

% generate input
input = zeros(1000, 1);
lib_util.randomgen(true);
for j=1:1:numel(input)
    input(j) = lib_util.randomgen() + 1i * lib_util.randomgen();
end

% run the resampler
output = lib_rx.resampling_polyphase(input, L, M, fir_coef_tx);
clear all;
close all;

% The newly introduced STF cover sequence requires some code changes.
% This script tests these changes when using a samplewise algorithm with accumulator + shift register.

% test a unch of different configurations
N_iter = 25;

% run multiple simulations
for iter = 1:N_iter

    % standard of random cover sequence
    cs_u1 = [1, -1, 1, 1, -1, -1, -1];
    cs_u248 = [1, -1, 1, 1, -1, -1, -1, -1, -1];
    cs_random_len = randi([4,41],1,1);
    cs_random = randi([0 1], cs_random_len, 1)*2-1;
    
    % pick one
    %cs = cs_u1;
    %cs = cs_u248;
    cs = cs_random;
    
    % random base sequence length with oversampling
    bs_len_os = randi([1,77],1,1);
    cs_os = repelem(cs, bs_len_os);
    
    % indexes of weight changes from 1 to -1 and -1 to 1 within the shift register due to the cover sequence
    acc_ptr_offset_weight_change_negative_vec = get_acc_ptr_offset_weight_change_negative(cs, bs_len_os);
    acc_ptr_offset_weight_change_positive_vec = get_acc_ptr_offset_weight_change_positive(cs, bs_len_os);
    
    % random IQ input
    inp_len = randi([101,9999],1,1);
    inp = randn(1, inp_len) + 1i*randn(1, inp_len);
    
    % reference algorithm convolution
    out_ref = conv(cs_os, inp);
    out_ref = out_ref(1:inp_len);
    
    % efficient (lowest number of multiplications, sequential memory access), samplewise (low latency) algorithm
    out_algo = zeros(size(out_ref));
    shift_register = zeros(size(cs_os));
    acc = 0;
    acc_ptr = 1;
    for i=1:inp_len
        % pop and push: adjust accumulator with values from the ends of the shift register
        acc = acc - cs(end)*shift_register(acc_ptr);
        acc = acc + cs(1)*inp(i);
    
        % adjust accumulator at weight changes within shift register
        for j=acc_ptr_offset_weight_change_negative_vec
            acc_ptr_offset = get_acc_ptr_relative(acc_ptr, numel(shift_register), j);
            acc = acc - 2*shift_register(acc_ptr_offset);
        end
        for j=acc_ptr_offset_weight_change_positive_vec
            acc_ptr_offset = get_acc_ptr_relative(acc_ptr, numel(shift_register), j);
            acc = acc + 2*shift_register(acc_ptr_offset);
        end
    
        % shift by overwriting value and ...
        shift_register(acc_ptr) = inp(i);
    
        % .. advancing and wrapping pointer
        acc_ptr = acc_ptr+1;
        if acc_ptr > numel(shift_register)
            acc_ptr = 1;
        end
        
        % save
        out_algo(i) = acc;
    end
    
    % compare
    real_err_max = max(abs(real(out_ref) - real(out_algo)));
    imag_err_max = max(abs(imag(out_ref) - imag(out_algo)));
    
    % show result
    %disp(real_err_max);
    %disp(imag_err_max);
    
    % warn user
    if real_err_max > 1e-8 || imag_err_max > 1e-8
        error("Not precise enough");
    end
end

function acc_ptr_offset_weight_change_negative_vec = get_acc_ptr_offset_weight_change_negative(cs, bs_len)
    
    acc_ptr_offset_weight_change_negative_vec = [];

    for i=2:numel(cs)
        if cs(i-1) == 1 && cs(i) == -1
            acc_ptr_offset_weight_change_negative_vec(end+1) = i-1;
        end
    end
    
    acc_ptr_offset_weight_change_negative_vec = acc_ptr_offset_weight_change_negative_vec*bs_len;
end

function acc_ptr_offset_weight_change_positive_vec = get_acc_ptr_offset_weight_change_positive(cs, bs_len)
    
    acc_ptr_offset_weight_change_positive_vec = [];

    for i=2:numel(cs)
        if cs(i-1) == -1 && cs(i) == 1
            acc_ptr_offset_weight_change_positive_vec(end+1) = i-1;
        end
    end
    
    acc_ptr_offset_weight_change_positive_vec = acc_ptr_offset_weight_change_positive_vec*bs_len;
end

function acc_ptr_offset = get_acc_ptr_relative(acc_ptr, acc_ptr_max, offset)

    acc_ptr_offset = acc_ptr - offset;

    if acc_ptr_offset <= 0
        acc_ptr_offset = acc_ptr_offset + acc_ptr_max;
    end
end

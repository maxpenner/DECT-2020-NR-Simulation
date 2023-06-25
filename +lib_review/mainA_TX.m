clear all;
close all;
clc;

% Script to compare TX implementation of Matlab and C++.
%
% Note: In Matlab, non-standard compliant features were inserted, which might have to be deactivated for a comporison against C++.
%
%   1) If symbols with STF show errors, go to file +lib_6_generic_procedures/ofdm_signal_generation_Cyclic_prefix_insertion.m and deactivate this line:
%
%       samples_antenna(1:STF_in_samples,:) = sqrt(N_eff_TX)*samples_antenna(1:STF_in_samples,:);
%
%   2) If symbols with DRS show errors, go to file +lib_5_physical_layer_transmissions/DRS.m and deactivate the following line:
%
%       y_b1 = y_b1*sqrt(N_eff_TX);
%
%   3) If symbols with PCC show errors, go to file +lib_6_generic_procedures/Transmit_diversity_precoding_Y.m and actiave this line
%
%       dont_use_preafactor_as_in_TS = true;
%
%      and deactivate this line:
%
%       dont_use_preafactor_as_in_TS = false;

warning on

% load and separate all json files
[filenames, n_files] = lib_util.get_all_filenames('results');
[tx_filenames, ~, ~] = lib_review.lib_helper.json_separate(filenames, n_files);

plot_every_packet = false;

% process each file
for i=1:1:numel(tx_filenames)

    % extract file folder and name
    file_struct = tx_filenames(i);
    filefolder = file_struct.folder;
    filename = file_struct.name;
    ffn = fullfile(filefolder, filename);
    disp(file_struct.name);

    % load json
    tx_json_struct = lib_review.lib_helper.json_load(ffn);

    % compare TX signals of C++ and Matlab numerically
    [~, ...
     samples_antenna_tx_matlab, ...
     samples_antenna_tx_matlab_resampled, ...
     samples_antenna_tx_cpp, ...
     dect_samp_rate, ...
     hw_samp_rate] = lib_review.lib_helper.TX_compare_numerically(tx_json_struct);

    % numerical comparison successful, now plot for visual confirmation
    if plot_every_packet == true
        lib_review.lib_helper.TX_compare_plot(  filename, ...
                                                samples_antenna_tx_matlab, ...
                                                samples_antenna_tx_matlab_resampled, ...
                                                samples_antenna_tx_cpp, ...
                                                dect_samp_rate, ...
                                                hw_samp_rate);
    end
end


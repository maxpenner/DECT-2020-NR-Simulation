clear all;
close all;
clc;

% This script is used to compare the C++ implementation.

% load all json filenames
[filenames, n_files] = lib_util.get_all_filenames('results');

% process each file
for i=1:1:n_files

    file_struct = filenames(i);

    disp(file_struct.name)

    % pick a test
    lib_review.packets_TX(file_struct.folder, file_struct.name);
    %lib_review.packets_CH(file_struct.folder, file_struct.name);
    %lib_review.packets_RX_FEC(file_struct.folder, file_struct.name);
end
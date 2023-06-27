clear all;
close all;
clc;

% script to compare TX and RX of C++ (exported with loopback firmware) with Matlab encoding/decoding

warning on

% load and separate all json files
[filenames, n_files] = lib_util.get_all_filenames('results');
[tx_filenames, ~, rx_synced_filenames] = lib_review.lib_helper.json_separate(filenames, n_files);

%% sync points from rx files
rx_synced_sync_fine_peak_time_64_vec = zeros(numel(rx_synced_filenames), 1);
for i=1:1:numel(rx_synced_filenames)

    % extract file folder and name
    file_struct = rx_synced_filenames(i);
    filefolder = file_struct.folder;
    filename = file_struct.name;
    ffn = fullfile(filefolder, filename);

    % load json
    rx_json_struct = lib_review.lib_helper.json_load(ffn);

    rx_synced_sync_fine_peak_time_64_vec(i) = rx_json_struct.meta.sync_fine_peak_time_64;
end

% find mapping between tx times and sync points
tx_to_sync_points_optimal_idx = lib_review.lib_helper.TX_compare_sync_points(   tx_filenames, ...
                                                                                "Sync points from rx files: ", ...
                                                                                rx_synced_sync_fine_peak_time_64_vec);

%% check every tx packet numerically, and then compare with rx packet

IQ_snr_dB_vec = ones(numel(tx_filenames), 1) * -1e9;
BER_PCC_vec = ones(numel(tx_filenames), 1) * -1e9;

plot_every_packet = false;

for i=1:1:numel(tx_filenames)

    if tx_to_sync_points_optimal_idx(i) <= 0
        continue;
    end

    % extract file folder and name
    file_struct = tx_filenames(i);
    filefolder = file_struct.folder;
    filename = file_struct.name;
    ffn = fullfile(filefolder, filename);
    disp(file_struct.name);

    % load tx json
    tx_json_struct = lib_review.lib_helper.json_load(ffn);

    % compare TX signals of C++ and Matlab numerically
    [tx, ...
     samples_antenna_tx_matlab, ...
     samples_antenna_tx_matlab_resampled, ...
     samples_antenna_tx_cpp, ...
     dect_samp_rate, ...
     hw_samp_rate] = lib_review.lib_helper.TX_compare_numerically(ffn, tx_json_struct);

    % extract RX file folder and name
    file_struct = rx_synced_filenames(tx_to_sync_points_optimal_idx(i));
    filefolder = file_struct.folder;
    filename = file_struct.name;
    ffn = fullfile(filefolder, filename);

    % load rx json
    rx_json_struct = lib_review.lib_helper.json_load(ffn);

    % first compare packets numerically
    [samples_antenna_rx_cpp, ...
     IQ_err_max_idx, ...
     IQ_snr_dB, ...
     BER_PCC] = lib_review.lib_helper.TX_RX_compare_numerically(tx, samples_antenna_tx_matlab, rx_json_struct);

    % then plot
    if plot_every_packet == true
        lib_review.lib_helper.TX_RX_compare_plot(samples_antenna_tx_matlab, samples_antenna_rx_cpp, IQ_err_max_idx);
    end

    % save for later plot
    IQ_snr_dB_vec(i) = IQ_snr_dB;
    BER_PCC_vec(i) = BER_PCC;
end

%% plot SNR

% remove any very low SNR
IQ_snr_dB_vec( IQ_snr_dB_vec < -1e-5) = [];

% plot and SNR histogramm
figure(2)
histogram(IQ_snr_dB_vec)
xlim([0 max(15, 1.5*max(IQ_snr_dB_vec))])
title_str = "SNR measured for " + num2str(numel(IQ_snr_dB_vec)) +  " files";
title(title_str);

%% plot BER of PCC

% remove any very low SNR
BER_PCC_vec( BER_PCC_vec < -1e-5) = [];

% plot and SNR histogramm
figure(3)
histogram(BER_PCC_vec,200)
xlim([-0.1 0.5])
xline(mean(BER_PCC_vec), 'r')
title_str = "BER of PCC " + num2str(numel(BER_PCC_vec)) +  " files";
title(title_str);

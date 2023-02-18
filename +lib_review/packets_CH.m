function [] = packets_CH(filefolder, filename)

%   It is fairly hard to test a DECT2020-NR receiver in C++ as is requires a proper channel model implementation.
%   This is particularly true for MIMO testcases with correlated, doubly selective channels.
%   Instead, this function can be used.   
%
%   1) With the help of the function lib_review.packets_TX, frames are generated and tested in C++.
%      Alternatively, random frames can be generate directly in Matlab.
%
%   2) The TX samples are passed through a channel simulation in Matlab.
%
%   3) The generated RX IQ are saved in a json, the IQ samples can then be read into C++ and processed there.
%

    verbose = 1;

    %% read samples with TX data, compare and return everything we need to create an RX signal
    [json_struct, mac_meta_tx, tx, samples_antenna_tx, resampling_param] = lib_review.packets_TX(filefolder, filename);
    
    %% create channel
    
    % number of antennas at TX and RX
    N_TX = tx.phy_4_5.tm_mode.N_TX;
    N_RX = 2;%size(samples_antenna_tx, 2);
    
    % RF channel parameters (see +lib_rf_channel/rf_channel.m)
    ch                      = lib_rf_channel.rf_channel();
    ch.verbose              = verbose;
    ch.verbose_cp           = tx.phy_4_5.numerology.N_b_CP*tx.mac_meta.oversampling;
    ch.type                 = 'deterministic';
    ch.amp                  = 1.0;
    ch.noise                = true;
    ch.snr_db               = 30;
    ch.spectrum_occupied    = tx.phy_4_5.n_spectrum_occupied/tx.mac_meta.oversampling;
    ch.N_TX                 = N_TX;
    ch.N_RX                 = N_RX;
    ch.awgn_random_source   = 'global';
    ch.awgn_randomstream    = RandStream('mt19937ar','Seed', randi(1e9,[1 1]));
    ch.d_sto                = 0;
    ch.d_cfo                = 0;
    ch.d_err_phase          = 0;
    ch.r_random_source      = 'global';
    ch.r_seed    	        = randi(1e9,[1 1]);
    ch.r_sto                = 0;
    ch.r_cfo                = 0;
    ch.r_err_phase          = 0;
    ch.r_samp_rate          = tx.phy_4_5.numerology.B_u_b_DFT*tx.mac_meta.oversampling;
    ch.r_max_doppler        = 1.946;                            % 1.946 19.458
    ch.r_type   	        = 'TDL-iii';
    ch.r_DS_desired         = 10^(-7.03 + 0.00*randn(1,1));
    ch.r_K                  = db2pow(9.0 + 0.00*randn(1,1));    %93e-9;
    ch.r_interpolation      = true;
    ch.r_gains_active 	    = true;
    %ch.init_rayleigh_rician_channel();
    
    %% pass tx signal through channel and yield rx signal (this step can be skipped)
    samples_antenna_rx = ch.pass_samples(samples_antenna_tx, 0);

    %% save data into json file
    
    % concatenate IQ samples
    samples_antenna_rx_vec = reshape(samples_antenna_rx, [], 1);

    % add RX IQ samples to json
    json_struct.dect.data.IQ.real = real(samples_antenna_rx_vec);
    json_struct.dect.data.IQ.imag = imag(samples_antenna_rx_vec);

    % receiver needs to know the correct number of antennas
    json_struct.dect.N_RX = N_RX;

    % save to file
    filename = 'rx_packet_' + convertCharsToStrings(sprintf( '%010d', json_struct.identifier));
    lib_review.save_json(fullfile('results/', filename), jsonencode(json_struct, PrettyPrint=true));
end


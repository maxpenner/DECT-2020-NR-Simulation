clear all;
close all;

% This script illustrates how to use this DECT-2020 New Radio link-level simulation environment.
% It defines a transmitter and a receiver handle class object.
% Also, it initiates the wireless channel model (this step can be skipped and replaced by a custom channel model).
% The transmitter generates a packet, sends it through the wireless channel and finally the receiver decodes it.

rng('shuffle');
%rng(1140598280);

warning('on');

% nice verbose plots
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

%profile on

%% these variables need to be set before creating tx
mac_meta_tx.u = 1;                                  % mu = 1, 2, 4 or 8
mac_meta_tx.b = 1;                                  % beta = 1, 2, 4, 8, 12 or 16
mac_meta_tx.PacketLengthType = 0;                   % 0 for subslots, 1 for slots
mac_meta_tx.PacketLength = 10;                      % min is 1, max is 16 according to Table 6.2.1-2a in part 4
mac_meta_tx.tm_mode_0_to_11 = 0;                    % Table 7.2-1, mode determines wether transmission is closed loop or not, values range from 0 to 11
mac_meta_tx.mcs_index = 1;                          % Table A-1 in part 3, values range from 0 to 11
mac_meta_tx.Z = 6144;                               % 5.3 -> so far only Z=6144 fully supported, 2048 only at TX, RX missing (Matlab has no option for Z=2048 in LTE toolbox)
mac_meta_tx.oversampling = 16;                      % By how much do we oversample our ofdm packet compared to critical sampling (insert zeros at specturm edges before IFFT)?
mac_meta_tx.codebook_index = 0;                     % 6.3.4, any value other than 0 makes packet beamformed, throws error if out of bound (depends on tm_mode_0_to_11)
mac_meta_tx.PLCF_type = 1;                          % Type 1 is 40 bits, Type 2 is 80 bits
mac_meta_tx.rv = 1;                                 % HARQ version, values range from 0, 1, 2 to 3 (right HARQ retransmission order is 0 2 3 1)
mac_meta_tx.network_id = de2bi(1e6,32,'left-msb');  % 7.6.6 must be given as a 32 bit vector with network_id(1) being the MSB, network_id must be known for scrambler on PHY

% temporary restrictions
if mac_meta_tx.Z ~= 6144
    error('Z must be 6144.');
end

%% create tx
verbose = 2;                                        % show data during execution: 0 false, 1 only text, 2 text + plots
tx = dect_tx(verbose, mac_meta_tx);

%% create rx

% assume the receiver has full knowledge of meta data at the transmitter (usually extracted from STF+PCC)
mac_meta_rx = mac_meta_tx;

% number of antennas at the receiver
mac_meta_rx.N_RX = 2;

% STO synchronization parameters (see +lib_rx/sync_STO.m)
mac_meta_rx.sto_config = lib_rx.sync_STO_param(mac_meta_tx.u, mac_meta_tx.b, mac_meta_tx.oversampling);
mac_meta_rx.sto_config.use_sto_sync = true;

% CFO synchronization parameters (see +lib_rx/sync_CFO.m)
mac_meta_rx.cfo_config = lib_rx.sync_CFO_param(mac_meta_tx.u);
mac_meta_rx.cfo_config.use_cfo_fractional = true;
mac_meta_rx.cfo_config.use_cfo_integer = true;
mac_meta_rx.cfo_config.use_cfo_residual = true;

% channel estimation parameters (see +lib_rx/channel_estimation_wiener.m)
mac_meta_rx.use_ch_estim_type = 'wiener';

% if we are using a wiener filter for channel esimation
if strcmp(mac_meta_rx.use_ch_estim_type,'wiener') == true

    % best case assumption for SNR, worst case assumptions for doppler and delay spread
    noise_estim_wiener               	= 1/10^(30/10);
    f_d_hertz_wiener                    = 100;
    tau_rms_sec_wiener                  = 300e-9;

    % the wiener filter coefficients
    wiener = lib_rx.channel_estimation_wiener_weights(  tx.phy_4_5.physical_resource_mapping_DRS_cell,...    % where are the pilots tones located?
                                                        tx.phy_4_5.numerology.N_b_DFT,...
                                                        tx.phy_4_5.N_PACKET_symb,...
                                                        tx.phy_4_5.numerology.N_b_CP, ...
                                                        tx.phy_4_5.numerology.B_u_b_DFT,...
                                                        noise_estim_wiener, ...
                                                        f_d_hertz_wiener, ...
                                                        tau_rms_sec_wiener);
end

% channel equalization (see dect_rx.m)
mac_meta_rx.use_equalization = true;

% create actual receiver
rx = dect_rx(verbose, mac_meta_rx);

% save receiver independent wiener coefficients
if strcmp(rx.mac_meta.use_ch_estim_type,'wiener') == true
    rx.wiener = wiener;
end

%% create channel

% number of antennas at TX and RX
N_TX = tx.phy_4_5.tm_mode.N_TX;
N_RX = rx.mac_meta.N_RX;

% RF channel parameters (see +lib_rf_channel/rf_channel.m)
ch                      = lib_rf_channel.rf_channel();
ch.verbose              = verbose;
ch.verbose_cp           = tx.phy_4_5.numerology.N_b_CP*tx.mac_meta.oversampling;
ch.type                 = 'rayleigh';
ch.amp                  = 1.0;
ch.noise                = true;
ch.snr_db               = 30;
ch.spectrum_occupied    = tx.phy_4_5.n_spectrum_occupied/tx.mac_meta.oversampling;
ch.N_TX                 = N_TX;
ch.N_RX                 = N_RX;
ch.awgn_random_source   = 'global';
ch.awgn_randomstream    = RandStream('mt19937ar','Seed', randi(1e9,[1 1]));
ch.d_sto                = 123 + 10*rx.mac_meta.sto_config.n_samples_STF_b_os;
ch.d_cfo                = 1.7*(1/(tx.phy_4_5.numerology.N_b_DFT*tx.mac_meta.oversampling));
ch.d_err_phase          = deg2rad(123);
ch.r_random_source      = 'global';
ch.r_seed    	        = randi(1e9,[1 1]);
ch.r_sto                = 123 + 10*rx.mac_meta.sto_config.n_samples_STF_b_os;
ch.r_cfo                = 13.75*(1/(tx.phy_4_5.numerology.N_b_DFT*tx.mac_meta.oversampling));
ch.r_err_phase          = deg2rad(123);
ch.r_samp_rate          = tx.phy_4_5.numerology.B_u_b_DFT*tx.mac_meta.oversampling;
ch.r_max_doppler        = 1.946;                            % 1.946 19.458
ch.r_type   	        = 'TDL-iii';
ch.r_DS_desired         = 10^(-7.03 + 0.00*randn(1,1));
ch.r_K                  = db2pow(9.0 + 0.00*randn(1,1));    %93e-9;
ch.r_interpolation      = true;
ch.r_gains_active 	    = true;
ch.init_rayleigh_rician_channel();

%% give rx handles so it can debug, e.g. perfect channel knowledge
rx.tx_handle = tx;
rx.ch_handle = ch;

%% generate tx signal

% generate random PCC bits
PCC_user_bits = [];
if mac_meta_tx.PLCF_type == 1
    PCC_user_bits = randi([0 1], 40, 1);
elseif mac_meta_tx.PLCF_type == 2
    PCC_user_bits = randi([0 1], 80, 1);
end

% how many PDC bits does tx need?
N_TB_bits = tx.phy_4_5.N_TB_bits;

% generate bits
PDC_user_bits = randi([0 1], N_TB_bits, 1);

% let tx create the packet
samples_antenna_tx = tx.generate_packet(PCC_user_bits, PDC_user_bits);

%% pass tx signal through channel and yield rx signal (this step can be skipped)
samples_antenna_rx = ch.pass_samples(samples_antenna_tx, 0);

%% let rx decode the frame
[PCC_user_bits_recovered, PDC_user_bits_recovered] = rx.demod_decode_packet(samples_antenna_rx);

%profile viewer
%profile off


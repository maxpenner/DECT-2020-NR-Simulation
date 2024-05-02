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
mac_meta_tx.PacketLength = 4;                       % min is 1, max is 16 according to Table 6.2.1-2a in part 4
mac_meta_tx.tm_mode_0_to_11 = 0;                    % Table 7.2-1, mode determines wether transmission is closed loop or not, values range from 0 to 11
mac_meta_tx.mcs_index = 5;                          % Table A-1 in part 3, values range from 0 to 11
mac_meta_tx.Z = 6144;                               % 5.3 -> so far only Z=6144 fully supported, 2048 only at TX, RX missing (Matlab has no option for Z=2048 in LTE toolbox)
mac_meta_tx.oversampling = 8;                       % By how much do we oversample our ofdm packet compared to critical sampling (insert zeros at specturm edges before IFFT)?
mac_meta_tx.codebook_index = 0;                     % 6.3.4, any value other than 0 makes packet beamformed, throws error if out of bound (depends on tm_mode_0_to_11)
mac_meta_tx.PLCF_type = 1;                          % Type 1 is 40 bits, Type 2 is 80 bits
mac_meta_tx.rv = 0;                                 % HARQ version, values range from 0, 1, 2 to 3 (right HARQ retransmission order is 0 2 3 1)
mac_meta_tx.network_id = de2bi(1e6,32,'left-msb');  % 7.6.6 must be given as a 32 bit vector with network_id(1) being the MSB, network_id must be known for scrambler on PHY

% temporary restrictions
if mac_meta_tx.Z ~= 6144
    error('Z must be 6144.');
end

%% create tx
verbose = 2;                                        % show data during execution: 0 false, 1 only text, 2 text + plots
tx = dect_tx(verbose, mac_meta_tx);

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

%% create rx

% assume the receiver has full knowledge of meta data at the transmitter (usually extracted from STF+PCC)
mac_meta_rx = mac_meta_tx;

% number of antennas at the receiver
mac_meta_rx.N_RX = 2;

% Synchronization based on STF (typically link-layer-simulations do not include synchronization)
%
% If synchronization is turned on (i.e. mac_meta_rx.synchronization.stf.active = true) the receiver class dect_rx will try
% synchronize a packet before decoding it. For that, the dect_rx member function demod_decode_packet(samples_antenna_rx) must
% be called with samples_antenna_rx having more samples than samples_antenna_tx.
%
% If synchronization is turned off (i.e. mac_meta_rx.synchronization.stf.active = false) the dect_rx member
% function demod_decode_packet(samples_antenna_rx) must be called with samples_antenna_tx having the exact same number
% of samples as samples_antenna_tx, but not necessarily the same number of antennas (depends on MIMO mode).
%
mac_meta_rx.synchronization.stf.active = true;
if mac_meta_rx.synchronization.stf.active == true

    % STO (detection, coarse peak search, fine peak search)
    mac_meta_rx.synchronization.stf.sto_config = lib_rx.sync_STO_param(mac_meta_tx.u, mac_meta_tx.b, mac_meta_tx.oversampling);
    
    % CFO (fractional, integer)
    mac_meta_rx.synchronization.stf.cfo_config = lib_rx.sync_CFO_param(mac_meta_tx.u);
    mac_meta_rx.synchronization.stf.cfo_config.active_fractional = true;
    mac_meta_rx.synchronization.stf.cfo_config.active_integer = true;
end

% synchronization of residual CFO based on DRS, not STF
mac_meta_rx.synchronization.drs.cfo_config.active_residual = true;

% create actual receiver
rx = dect_rx(verbose, mac_meta_rx);

%% create channel, can be replaced with a custom channel type

% number of antennas at TX and RX
N_TX = tx.phy_4_5.tm_mode.N_TX;
N_RX = rx.mac_meta.N_RX;

% RF channel parameters (see +lib_rf_channel/rf_channel.m) valid for all channel types.
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

% Parameters with d_ only used if ch.type = 'deterministic'. It is called deterministic because STO, CFO and error phase have fixed, known values.
ch.d_sto                = 123 + 2*size(samples_antenna_tx, 1);
ch.d_cfo                = 1.7*(1/(tx.phy_4_5.numerology.N_b_DFT*tx.mac_meta.oversampling));
ch.d_err_phase          = deg2rad(123);

% Parameters with r_ only used if ch.type = 'rayleigh' or 'rician'.
ch.r_random_source      = 'global';
ch.r_seed    	        = randi(1e9,[1 1]);
ch.r_sto                = ch.d_sto;
ch.r_cfo                = ch.d_cfo;
ch.r_err_phase          = ch.d_err_phase;
ch.r_samp_rate          = tx.phy_4_5.numerology.B_u_b_DFT*tx.mac_meta.oversampling;
ch.r_max_doppler        = 1.946;                            % 1.946 19.458
ch.r_type   	        = 'TDL-iii';
ch.r_DS_desired         = 10^(-7.03 + 0.00*randn(1,1));
ch.r_K                  = db2pow(9.0 + 0.00*randn(1,1));    %93e-9;
ch.r_interpolation      = true;
ch.r_gains_active       = true;
ch.init_rayleigh_rician_channel();

%% give rx handles so it can debug, e.g. perfect channel knowledge
rx.tx_handle = tx;
rx.ch_handle = ch;

%% pass tx signal through channel and yield rx signal (this step can be skipped)
samples_antenna_rx = ch.pass_samples(samples_antenna_tx, 0);

%% The property snr_db of the RF channel refers to the inband noise. If oversampling is used, we have to remove out-of-band noise, otherwise synchronization in time domain is impaired.
if mac_meta_tx.oversampling > 1
    % Kaiser LPF
    lowpassfilter = designfilt( 'lowpassfir', ...
                                'PassbandFrequency', 0.6 / mac_meta_tx.oversampling, ...
                                'StopbandFrequency',0.8 / mac_meta_tx.oversampling, ...
                                'PassbandRipple', 10, ...
                                'StopbandAttenuation', 30, ...
                                'SampleRate', 1, ...
                                'DesignMethod', 'kaiserwin', ...
                                'MinOrder', 'even');
    
    assert(mod(numel(lowpassfilter.Coefficients), 2) == 1, 'fractional LPF delay');
    
    samples_antenna_rx = filter(lowpassfilter, samples_antenna_rx);
    
    % compensate for deterministic filter delay prior to synchronization
    lowpassfilter_delay = (numel(lowpassfilter.Coefficients)-1)/2;
    samples_antenna_rx(1:end-lowpassfilter_delay) = samples_antenna_rx(lowpassfilter_delay+1:end);
end

%% let rx decode the frame
[PCC_user_bits_recovered, PDC_user_bits_recovered] = rx.demod_decode_packet(samples_antenna_rx);

%profile viewer
%profile off

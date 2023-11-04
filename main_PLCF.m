clear all;
close all;
clc;

% This script generates the Physical Header Field of DECT-2020 New Radio packets.
% The length is 40 or 80 bits.
% Refer to Table 6.2.1-1, 6.2.1-2 and 6.2.1-2a of part 4 of the TS.

%% these variables need to be set before creating tx
mac_meta_tx.u = 1;                                  % mu = 1, 2, 4 or 8
mac_meta_tx.b = 1;                                  % beta = 1, 2, 4, 8, 12 or 16
mac_meta_tx.PacketLengthType = 0;                   % 0 for subslots, 1 for slots
mac_meta_tx.PacketLength = 2;                       % min is 1, max is 16 according to Table 6.2.1-2a in part 4
mac_meta_tx.tm_mode_0_to_11 = 0;                 	% Table 7.2-1, mode determines wether transmission is closed loop or not, values range from 0 to 11
mac_meta_tx.mcs_index = 1;                          % Table A-1 in part 3, values range from 0 to 11
mac_meta_tx.Z = 6144;                               % 5.3 -> so far only Z=6144 fully supported, 2048 only at TX, RX missing (Matlab has no option for Z=2048 in LTE toolbox)
mac_meta_tx.oversampling = 16;                    	% By how much do we oversample our ofdm packet compared to critical sampling (insert zeros at specturm edges before IFFT)?
mac_meta_tx.codebook_index = 0;                     % 6.3.4, any value other than 0 makes packet beamformed, throws error if out of bound (depends on tm_mode_0_to_11)
mac_meta_tx.PLCF_type = 2;                          % Type 1 is 40 bits, Type 2 is 80 bits
mac_meta_tx.rv = 1;                                 % HARQ version, values range from 0, 1, 2 to 3 (right HARQ retransmission order is 0 2 3 1)
mac_meta_tx.network_id = de2bi(1e6,32,'left-msb');  % 7.6.6 must be given as a 32 bit vector with network_id(1) being the MSB, network_id must be known for scrambler on PHY

% temporary restrictions
if mac_meta_tx.Z ~= 6144
    error('Z must be 6144.');
end

%% Generate real PCC based on meta data, errors are thrown for any value out of bound.

% Refer to Table 6.2.1-1, 6.2.1-2 and 6.2.1-2a of part 4 of the TS.

plcf_meta.type                  = mac_meta_tx.PLCF_type;
plcf_meta.Header_Format         = 0;
plcf_meta.Packet_length_type    = mac_meta_tx.PacketLengthType;
plcf_meta.Packet_length         = mac_meta_tx.PacketLength - 1;     % Table 6.2.1-1: "signalled numerical value plus one"
plcf_meta.Short_Network_ID      = 123;
plcf_meta.Transmitter_Identity  = 123;
plcf_meta.Transmit_Power        = 0;
plcf_meta.DF_MCS                = mac_meta_tx.mcs_index + 1;        % Table 6.2.2-3, effectively signalled numerical value plus one

if plcf_meta.type == 2

    plcf_meta.Receiver_Identity = 123;

    switch mac_meta_tx.tm_mode_0_to_11
        case {0,1,3,5,7,10}
            plcf_meta.Number_of_Spatial_Streams = 0; % do not change
        case {2,4,8}
            plcf_meta.Number_of_Spatial_Streams = 1; % do not change
        case {6,9}
            plcf_meta.Number_of_Spatial_Streams = 2; % do not change
        case 11
            plcf_meta.Number_of_Spatial_Streams = 3; % do not change
        otherwise
            error('Unknown transmission type');
    end

    % only Table 6.2.1-2
    if plcf_meta.Header_Format == 0
        plcf_meta.DF_Redundancy_Version      = mac_meta_tx.rv;
        plcf_meta.DF_New_data_Indication     = 0;
        plcf_meta.DF_HARQ_Process_Number     = 0;
    end

    plcf_meta.Feedback_format = 1;
    switch plcf_meta.Feedback_format
        case 1
            plcf_meta.HARQ_Process_number    = 0;
            plcf_meta.Transmission_feedback  = 0;
            plcf_meta.CQI                    = 0;
            plcf_meta.Buffer_Status          = 0;
        case 2
            plcf_meta.CQI                    = 0;
            plcf_meta.Buffer_Status          = 0;
            plcf_meta.MIMO_Feedback          = 0;
            plcf_meta.Codebook_Index         = 0;
        case 3
            plcf_meta.HARQ_Process_number0   = 0;
            plcf_meta.Transmission_feedback0 = 0;
            plcf_meta.HARQ_Process_number1   = 0;
            plcf_meta.Transmission_feedback1 = 0;
            plcf_meta.CQI                    = 0;
        case 4
            plcf_meta.HARQ_feedback_bitmap   = 0;
            plcf_meta.CQI                    = 0;
        case 5
            plcf_meta.HARQ_Process_number    = 0;
            plcf_meta.Transmission_feedback  = 0;
            plcf_meta.MIMO_Feedback          = 0;
            plcf_meta.Codebook_Index         = 0;
        otherwise
            error('Feedback format must be 1, 2, 3, 4 or 5.');
    end
end

PCC_user_bits = lib_part4.Physical_Header_Field(plcf_meta);


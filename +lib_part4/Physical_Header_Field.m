function [PCC_user_bits] = Physical_Header_Field(plcf_meta)

    % Definition of plcf.
    % Refer to Table 6.2.1-1, 6.2.1-2 and 6.2.1-2a of part 4 of the TS.

    plcf.type = plcf_meta.type;
    
    % the length is defined
    if plcf.type == 1
        plcf.length = 40;
    elseif plcf.type == 2
        plcf.length = 80;
    else
        error('Type must be 1 or 2.');
    end

    % next we have to choose a header format
    if plcf.type == 1

        if plcf_meta.Header_Format ~= 0
            error('Header format must be 0 for type 1.');
        end

        plcf.Header_Format = de2bi(0, 3,'left-msb');

    elseif plcf.type == 2

        if plcf_meta.Header_Format == 0
            plcf.Header_Format = de2bi(0, 3,'left-msb');    % must be chosen when using HARQ
        elseif plcf_meta.Header_Format == 1
            plcf.Header_Format = de2bi(1, 3,'left-msb');
        else
            error('Header format must be 0 or 1 for type 2.');
        end
    end
    
    % the other data
    % all tables
    plcf.Packet_length_type     = de2bi(plcf_meta.Packet_length_type,   1,'left-msb');
    plcf.Packet_length          = de2bi(plcf_meta.Packet_length,        4,'left-msb');
    plcf.Short_Network_ID       = de2bi(plcf_meta.Short_Network_ID,     8,'left-msb');
    plcf.Transmitter_Identity   = de2bi(plcf_meta.Transmitter_Identity,16,'left-msb');
    plcf.Transmit_Power         = de2bi(plcf_meta.Transmit_Power,       4,'left-msb');
    
    % only Table 6.2.1-1
    if plcf.type == 1
        
        plcf.Reserved                   = de2bi(0,1,'left-msb');
        plcf.DF_MCS                     = de2bi(plcf_meta.DF_MCS,3,'left-msb');
        
    % only Table 6.2.1-2 and Table 6.2.1-2a
    elseif plcf.type == 2

        plcf.DF_MCS                     = de2bi(plcf_meta.DF_MCS,                    4,'left-msb');	% Table 6.2.2-3
        plcf.Receiver_Identity          = de2bi(plcf_meta.Receiver_Identity,        16,'left-msb');
        plcf.Number_of_Spatial_Streams  = de2bi(plcf_meta.Number_of_Spatial_Streams, 2,'left-msb');

        % only Table 6.2.1-2 
        if sum(plcf.Header_Format) == 0
            
            plcf.DF_Redundancy_Version      = de2bi(plcf_meta.DF_Redundancy_Version,    2,'left-msb');
            plcf.DF_New_data_Indication     = de2bi(plcf_meta.DF_New_data_Indication,   1,'left-msb');
            plcf.DF_HARQ_Process_Number     = de2bi(plcf_meta.DF_HARQ_Process_Number,   3,'left-msb');
            
        % only Table 6.2.1-2a
        elseif sum(plcf.Header_Format) == 1
            plcf.Reserved                   = de2bi(0,6,'left-msb');
        else
            error('Incorrect header format.');
        end

        plcf.Feedback_format                = de2bi(plcf_meta.Feedback_format,      4,'left-msb');    % Table 6.2.2-1
        switch bi2de(plcf.Feedback_format,'left-msb')
            case 1
                plcf.HARQ_Process_number    = de2bi(plcf_meta.HARQ_Process_number,  3,'left-msb');
                plcf.Transmission_feedback  = de2bi(plcf_meta.Transmission_feedback,1,'left-msb');
                plcf.CQI                    = de2bi(plcf_meta.CQI,                  4,'left-msb');
                plcf.Buffer_Status          = de2bi(plcf_meta.Buffer_Status,        4,'left-msb');
                plcf.Feedback_info          = [plcf.HARQ_Process_number, plcf.Transmission_feedback, plcf.CQI, plcf.Buffer_Status];
            case 2
                plcf.CQI                    = de2bi(plcf_meta.CQI,                  4,'left-msb');
                plcf.Buffer_Status          = de2bi(plcf_meta.Buffer_Status,        4,'left-msb');
                plcf.MIMO_Feedback          = de2bi(plcf_meta.MIMO_Feedback,        1,'left-msb');
                plcf.Codebook_Index         = de2bi(plcf_meta.Codebook_Index,       3,'left-msb');
                plcf.Feedback_info          = [plcf.CQI, plcf.Buffer_Status, plcf.MIMO_Feedback, plcf.Codebook_Index];
            case 3
                plcf.HARQ_Process_number0   = de2bi(plcf_meta.HARQ_Process_number0,      3,'left-msb');
                plcf.Transmission_feedback0 = de2bi(plcf_meta.Transmission_feedback0,    1,'left-msb');
                plcf.HARQ_Process_number1   = de2bi(plcf_meta.HARQ_Process_number1,      3,'left-msb');
                plcf.Transmission_feedback1 = de2bi(plcf_meta.Transmission_feedback1,    1,'left-msb');
                plcf.CQI                    = de2bi(plcf_meta.CQI,                       4,'left-msb');
                plcf.Feedback_info          = [ plcf.HARQ_Process_number0, plcf.Transmission_feedback0,...
                                                plcf.HARQ_Process_number1, plcf.Transmission_feedback1,...
                                                plcf.CQI];
            case 4
                plcf.HARQ_feedback_bitmap   = de2bi(plcf_meta.HARQ_feedback_bitmap,      8,'left-msb');     
                plcf.CQI                    = de2bi(plcf_meta.CQI,                       4,'left-msb');
                plcf.Feedback_info          = [plcf.HARQ_feedback_bitmap, plcf.CQI];            
            case 5
                plcf.HARQ_Process_number    = de2bi(plcf_meta.HARQ_Process_number,       3,'left-msb');
                plcf.Transmission_feedback  = de2bi(plcf_meta.Transmission_feedback,     1,'left-msb');
                plcf.MIMO_Feedback          = de2bi(plcf_meta.MIMO_Feedback,             2,'left-msb');
                plcf.Codebook_Index         = de2bi(plcf_meta.Codebook_Index,            6,'left-msb');
                plcf.Feedback_info          = [plcf.HARQ_Process_number, plcf.Transmission_feedback, plcf.MIMO_Feedback, plcf.Codebook_Index];            
            otherwise
                error('Incorrect Feedback info format.');
        end
    else
        error('Incorrect plcf type.');
    end

    % generate actual PCC bits
    if plcf.type == 1
        
        PCC_user_bits = [   plcf.Header_Format,...
                            plcf.Packet_length_type,...
                            plcf.Packet_length,...
                            plcf.Short_Network_ID,...
                            plcf.Transmitter_Identity,...
                            plcf.Transmit_Power,...
                            plcf.Reserved,...
                            plcf.DF_MCS];
                        
        % sanity check
        if numel(PCC_user_bits) ~= 40
            error('Incorrect number of bits in plcf a.k.a. PCC user bits.');
        end
        
    % only Table 6.2.1-2 and Table 6.2.1-2a
    elseif plcf.type == 2
        
        if sum(plcf.Header_Format) == 0
            
            PCC_user_bits = [   plcf.Header_Format,...
                                plcf.Packet_length_type,...
                                plcf.Packet_length,...
                                plcf.Short_Network_ID,...
                                plcf.Transmitter_Identity,...
                                plcf.Transmit_Power,...
                                plcf.DF_MCS,...
                                plcf.Receiver_Identity,...
                                plcf.Number_of_Spatial_Streams,...
                                plcf.DF_Redundancy_Version,...
                                plcf.DF_New_data_Indication,...
                                plcf.DF_HARQ_Process_Number,...
                                plcf.Feedback_format,...
                                plcf.Feedback_info];
                        
            % sanity check
            if numel(PCC_user_bits) ~= 80
                error('Incorrect number of bits in plcf a.k.a. PCC user bits.');
            end
            
        % only Table 6.2.1-2a
        elseif sum(plcf.Header_Format) == 1
            
            PCC_user_bits = [   plcf.Header_Format,...
                                plcf.Packet_length_type,...
                                plcf.Packet_length,...
                                plcf.Short_Network_ID,...
                                plcf.Transmitter_Identity,...
                                plcf.Transmit_Power,...
                                plcf.DF_MCS,...
                                plcf.Receiver_Identity,...
                                plcf.Number_of_Spatial_Streams,...
                                plcf.Reserved,...
                                plcf.Feedback_format,...
                                plcf.Feedback_info];
                        
            % sanity check
            if numel(PCC_user_bits) ~= 80
                error('Incorrect number of bits in plcf a.k.a. PCC user bits.');
            end
        end
    end

    % sanity check
    if numel(PCC_user_bits) ~= plcf.length
        error('Incorrect PLCF size.');
    end
end


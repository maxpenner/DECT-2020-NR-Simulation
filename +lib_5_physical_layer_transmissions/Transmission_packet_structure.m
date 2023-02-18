function [N_PACKET_symb] = Transmission_packet_structure(numerology, PacketLengthType, PacketLength, N_eff_TX)

    switch PacketLengthType
        case 0
            N_PACKET_symb = PacketLength * numerology.N_SLOT_u_symb / numerology.N_SLOT_u_subslot;
        case 1
            N_PACKET_symb = PacketLength * numerology.N_SLOT_u_symb;
        otherwise
            error('PacketLengthType is %d, but it must be 0 or 1.', PacketLengthType);
    end

    if N_eff_TX >= 4

%         % according to standard we need only 15 symbols
%         if N_PACKET_symb < 15
%             error('For N_eff_TX >= 4 at least 15 symbols are required, see 5.1 Transmission packet structure.');
%         end

        % this is probably an error in the TS, we actually need at least 15 or 20 OFDM symbols
        if N_eff_TX == 4
            if N_PACKET_symb < 15
                error('For N_eff_TX == 4 at least 15 symbols are required, see 5.1 Transmission packet structure.');
            end
        elseif N_eff_TX == 8
            if N_PACKET_symb < 20 || mod(N_PACKET_symb, 10) ~= 0
                error('For N_eff_TX == 8 at least 20 symbols are required and N_PACKET_symb must be a multiple of 10, see 5.1 Transmission packet structure.');
            end
        else
            error("N_eff_TX unknown");
        end
    end
    
    % NOTE
    % ?

    % Figures 4.5-2 and 4.5-3
    tx_packet_len = N_PACKET_symb*numerology.T_u_symb;
    tx_packet_len_STF = 14/9*numerology.T_u_symb;
    tx_packet_len_DF = (N_PACKET_symb-2)*numerology.T_u_symb;
    tx_packet_len_GI = 4/9*numerology.T_u_symb;

    % sanity check
    %assert(tx_packet_len == tx_packet_len_STF+tx_packet_len_DF+tx_packet_len_GI);
    lib_dbg.check_equality_of_two_numbers(tx_packet_len, tx_packet_len_STF+tx_packet_len_DF+tx_packet_len_GI, 1e-6);
end


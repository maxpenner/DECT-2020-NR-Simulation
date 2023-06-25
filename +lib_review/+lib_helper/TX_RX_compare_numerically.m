function [samples_antenna_rx_cpp, ...
          IQ_err_max_idx, ...
          IQ_snr_dB, ...
          BER_PCC] ...
              = TX_RX_compare_numerically(  tx, ...
                                            samples_antenna_tx_matlab_resampled, ...
                                            rx_json_struct)

    N_RX = rx_json_struct.N_RX;

    %% compare IQ samples

    % read tx samples from json file from C++
    samples_antenna_rx_cpp = rx_json_struct.data.IQ.real + 1i*rx_json_struct.data.IQ.imag;

    % json concatenates antenna samples, reshape into matrix of size N_TX x N_samples_per_antenna
    N_samples_per_antenna = numel(samples_antenna_rx_cpp)/N_RX;
    samples_antenna_rx_cpp = reshape(samples_antenna_rx_cpp, N_samples_per_antenna, N_RX);

    % find smaller of both sizes
    [N_samples_matlab, N_TX_matlab] = size(samples_antenna_tx_matlab_resampled);
    [N_samples_cpp, N_RX_matlab] = size(samples_antenna_rx_cpp);
    N_samples_min = min(N_samples_matlab, N_samples_cpp);

    % this comparison only works for SISO
    if N_TX_matlab ~= 1 || N_TX_matlab ~= N_RX_matlab
        error("Comparison of TX signal with RX signal only possible for SISO.");
    end

    % find largest error between samples
    IQ_err = samples_antenna_tx_matlab_resampled(1:N_samples_min) - samples_antenna_rx_cpp(1:N_samples_min);
    [IQ_err_max, IQ_err_max_idx] = max(abs(IQ_err));
    fprintf("\tMaximum deviation between TX signal and RX signal is %f and idx=%d.\n", IQ_err_max, IQ_err_max_idx);

    if IQ_err_max > 1000.0
        error("TX samples deviate too much from RX samples.");
    end

    %% calculate SNR based on IQ samples
    err_power = rms(IQ_err)^2;
    IQ_snr_dB = 10*log10(1/err_power);
    fprintf("\tSNR between TX signal and RX signal is %f dB.\n", IQ_snr_dB);

    %% compare PCC

    % first extract bits from rx
    pcc_rx_as_short = rx_json_struct.data.binary.PCC;

    % convert to +1 and -1
    pcc_rx = sign(pcc_rx_as_short);
    pcc_tx = 2*double(tx.packet_data.pcc_enc_dbg.d) - 1;

    err_bits_pcc = abs(pcc_rx - pcc_tx);
    BER_PCC = sum(err_bits_pcc) / numel(err_bits_pcc);

    fprintf("\tBER of PCC is %f.\n", BER_PCC);
end
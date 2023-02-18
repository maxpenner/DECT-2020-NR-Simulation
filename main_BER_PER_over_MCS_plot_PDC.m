clear all;
close all;

% This scripts plots the results generated with main_BER_PER_over_MCS.m, main_BER_PER_over_MCS_convergence.m or main_BER_PER_over_MCS_convergence_full.m.

load('results/var_all.mat');

ber_global = n_bits_PDC_error_global./n_bits_PDC_sent_global;
per_global = n_packets_PDC_error_global./n_packets_PDC_sent_global;

% required plot configuration
K = db2pow(9.0);
use_export_fig_png = false;
use_export_fig_eps = false;

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

% plot configuration
colors = [0, 0.4470, 0.7410;...
            0.8500, 0.3250, 0.0980;...
            0.9290, 0.6940, 0.1250;...
            0.4940, 0.1840, 0.5560;...
            0.4660, 0.6740, 0.1880;...
            0.3010, 0.7450, 0.9330;...
            0.6350, 0.0780, 0.1840;...
            0.6350, 0.0780, 0.1840;...
            0.6350, 0.0780, 0.1840;...
            0.6350, 0.0780, 0.1840;...
            0.6350, 0.0780, 0.1840;...
            0.6350, 0.0780, 0.1840;...
            0.6350, 0.0780, 0.1840;...
            0.6350, 0.0780, 0.1840];
legend_font_size = 8;   
marker_size = 4;
axis_lim = [-15 55 1e-7 1e1];
legend_location = 'NorthEast';

figure()
clf()
for cnt = 1:1:numel(mcs_index_vec)

    % bits per symbol
    bps = bps_global(cnt);
    M = 2^bps;

    % convert subcarrier snr to ebn0
    %
    %   S/N = R_b * E_b / (B * N0)
    %
    %   R_b = bps / T
    %   B = 1 / T
    %
    %   S/N = bps * E_b/N0
    %
    % source: https://www.gaussianwaves.com/2008/11/relation-between-ebn0-and-snr-2/
    EbN0_vec = db2pow(snr_db_vec_global(cnt,:))/bps;
    EbN0_dB_vec = pow2db(EbN0_vec);

        % AWGN, diversity order 1
        if M==2
            ber_rayleigh = berawgn(EbN0_dB_vec,'psk',M,'nondiff');
        else
            ber_rayleigh = berawgn(EbN0_dB_vec,'qam',M);
        end
        str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)), ', AWGN');
        semilogy(snr_db_vec_global(cnt,:), ber_rayleigh,'-','DisplayName',str, 'Color', colors(cnt,:));
        hold on

            % RAYLEIGH, diversity order 1
            if M==2
                ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,1);
            else
                ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,1);
            end
            str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)), ', rayleigh');
            semilogy(snr_db_vec_global(cnt,:), ber_rayleigh,'--','DisplayName',str, 'Color', colors(cnt,:));
            hold on

                % RAYLEIGH, diversity order 2
                if M==2
                    ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,2);
                else
                    ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,2);
                end
                str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)), ', rayleigh div=2');
                semilogy(snr_db_vec_global(cnt,:), ber_rayleigh,'-.','DisplayName',str, 'Color', colors(cnt,:));
                hold on

                    % RICIAN, diversity order 1
                    if M==2
                        ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,1,K);
                    else
                        ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,1,K);
                    end
                    str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)), ', rician');
                    semilogy(snr_db_vec_global(cnt,:), ber_rayleigh,'-o','DisplayName',str, 'Color', colors(cnt,:));
                    hold on

                        % RICIAN, diversity order 1
                        if M==2
                            ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,2,K);
                        else
                            ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,2,K);
                        end
                        str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)), ', rician div=2');
                        semilogy(snr_db_vec_global(cnt,:), ber_rayleigh,'-d','DisplayName',str, 'Color', colors(cnt,:));
                        hold on

    % SIMULATION

    % our measured ber
    str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)));
    lineH = semilogy(snr_db_vec_global(cnt,:), ber_global(cnt,:), '-.o','DisplayName',str, 'Color', colors(cnt,:), 'MarkerSize', marker_size, 'MarkerFaceColor', colors(cnt,:));
    hold on

    % general settings
    title('BER uncoded')
    xlabel('SNR (dB)')
    ylabel('BER uncoded')
    legend('Location',legend_location, 'FontSize', legend_font_size)
    grid on
    axis(axis_lim)
    set(gca, 'ColorOrder', jet(100))
end

savefig('results/D_BER_SNR_PDC.fig')
if use_export_fig_png == true
    addpath('export_fig_lib')
    export_fig results/D_BER_SNR_PDC.png -m4 -transparent
end
if use_export_fig_eps == true
    addpath('export_fig_lib')
    export_fig results/D_BER_SNR_PDC.eps -m2.5 -transparent
end

figure()
clf()
for cnt = 1:1:numel(mcs_index_vec)

    % bits per symbol
    bps = bps_global(cnt);
    M = 2^bps;

    % convert subcarrier snr to ebn0
    %
    %   S/N = R_b * E_b / (B * N0)
    %
    %   R_b = bps / T
    %   B = 1 / T
    %
    %   S/N = bps * E_b/N0
    %
    % source: https://www.gaussianwaves.com/2008/11/relation-between-ebn0-and-snr-2/
    EbN0_vec = db2pow(snr_db_vec_global(cnt,:))/bps;
    EbN0_dB_vec = pow2db(EbN0_vec);

        % AWGN, diversity order 1
        if M==2
            ber_rayleigh = berawgn(EbN0_dB_vec,'psk',M,'nondiff');
        else
            ber_rayleigh = berawgn(EbN0_dB_vec,'qam',M);
        end
        str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)), ', AWGN');
        semilogy(EbN0_dB_vec, ber_rayleigh,'-','DisplayName',str, 'Color', colors(cnt,:));
        hold on

            % RAYLEIGH, diversity order 1
            if M==2
                ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,1);
            else
                ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,1);
            end
            str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)), ', rayleigh');
            semilogy(EbN0_dB_vec, ber_rayleigh,'--','DisplayName',str, 'Color', colors(cnt,:));
            hold on

                % RAYLEIGH, diversity order 2
                if M==2
                    ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,2);
                else
                    ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,2);
                end
                str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)), ', rayleigh div=2');
                semilogy(EbN0_dB_vec, ber_rayleigh,'-.','DisplayName',str, 'Color', colors(cnt,:));
                hold on

                    % RICIAN, diversity order 1
                    if M==2
                        ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,1,K);
                    else
                        ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,1,K);
                    end
                    str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)), ', rician');
                    semilogy(EbN0_dB_vec, ber_rayleigh,'-o','DisplayName',str, 'Color', colors(cnt,:));
                    hold on

                        % RICIAN, diversity order 1
                        if M==2
                            ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,2,K);
                        else
                            ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,2,K);
                        end
                        str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)), ', rician div=2');
                        semilogy(EbN0_dB_vec, ber_rayleigh,'-d','DisplayName',str, 'Color', colors(cnt,:));
                        hold on

    % SIMULATION

    % our measured ber
    str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)));
    lineH = semilogy(EbN0_dB_vec, ber_global(cnt,:), '-.o','DisplayName',str, 'Color', colors(cnt,:), 'MarkerSize', marker_size, 'MarkerFaceColor', colors(cnt,:));
    hold on

    % general settings
    title('BER uncoded')
    xlabel('EbN0 (dB)')
    ylabel('BER uncoded')
    legend('Location',legend_location, 'FontSize', legend_font_size)
    grid on
    axis(axis_lim)
    set(gca, 'ColorOrder', jet(100))
end

savefig('results/E_BER_EBN0_PDC.fig')
if use_export_fig_png == true
    addpath('export_fig_lib')
    export_fig results/E_BER_EBN0_PDC.png -m4 -transparent
end
if use_export_fig_eps == true
    addpath('export_fig_lib')
    export_fig results/E_BER_EBN0_PDC.eps -m2.5 -transparent
end

figure()
clf()
for cnt = 1:1:numel(mcs_index_vec)
    str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)));
    semilogy(snr_db_vec_global(cnt,:), per_global(cnt,:),'-o','DisplayName',str, 'Color', colors(cnt,:), 'MarkerSize', marker_size, 'MarkerFaceColor', colors(cnt,:));
    hold on
    title('Packet Decoding Performance')
    xlabel('SNR (dB)')
    ylabel('PER')
    legend('Location',legend_location, 'FontSize', legend_font_size)
    grid on
    axis(axis_lim)
end

savefig('results/F_PER_SNR_PDC.fig')
if use_export_fig_png == true
    addpath('export_fig_lib')
    export_fig results/F_PER_SNR_PDC.png -m4 -transparent
end
if use_export_fig_eps == true
    addpath('export_fig_lib')
    export_fig results/F_PER_SNR_PDC.eps -m2.5 -transparent
end

% ####################
% ####################
% save file with values
% ####################
% ####################

% replace NaN with zero
per_global(isnan(per_global)) = 0;

if exist('results/PER_SNR_PDC.txt', 'file') == 2
  delete('results/PER_SNR_PDC.txt');
end
clc
diary results/PER_SNR_PDC.txt
fprintf('%1.8f\n',snr_db_vec_global);
disp(' ');
disp(' ');
fprintf('%1.8f\n',per_global);
disp(' ');
disp(' ');
fprintf('%1.8f, ',snr_db_vec_global);
disp(' ');
disp(' ');
fprintf('%1.8f, ',per_global);
disp(' ');
disp(' ');
diary off

%% STF synchronization plots
if exist('STO_hist_coarse_global','var') == true

    % how many packets did we create per SNR?
    n_packets_per_SNR_coarse = sum(STO_hist_coarse_global,1);
    n_packets_per_SNR_fine = sum(STO_hist_fine_global,1);

    % sanity check
    if sum(n_packets_per_SNR_coarse - n_packets_per_SNR_fine) > 0
        error('Number of packets no equal.');
    end

    % how many packets are outside the evaluation range (outer edges)?
    n_packets_per_missed_coarse = STO_hist_coarse_global(1,:) + STO_hist_coarse_global(end,:);
    n_packets_per_missed_fine = STO_hist_fine_global(1,:) + STO_hist_fine_global(end,:);

    % we need a percentage of missed packets
    perc_packets_missed_coarse = n_packets_per_missed_coarse./n_packets_per_SNR_coarse;
    perc_packets_missed_fine = n_packets_per_missed_fine./n_packets_per_SNR_fine;

    % size of plot
    [n_bins, n_SNR_steps] = size(STO_hist_coarse_global);

    figure();
    clf()

    for i=1:1:n_SNR_steps

        row = (i-1)*2+1;

        % COARSE

        subplot(n_SNR_steps, 2, row);
        histogram('BinEdges',edges,'BinCounts',STO_hist_coarse_global(:,i), 'Normalization','pdf');

        xlim([edges(2) edges(end-1)]);
        ylim([0 1.1]);

        xlabel('STO Error in oversampled Samples')

        title_string = "SNR=" + num2str(snr_db_vec_global(i)) + "dB   Coarse";
        title(title_string);

        text(0,0.8,"Missed=" + num2str(perc_packets_missed_coarse(i)));

        % FINE

        subplot(n_SNR_steps, 2, row+1);
        histogram('BinEdges',edges,'BinCounts',STO_hist_fine_global(:,i), 'Normalization','pdf');

        xlim([edges(2) edges(end-1)]);
        ylim([0 1.1]);

        xlabel('STO Error in oversampled Samples')

        title_string = "SNR=" + num2str(snr_db_vec_global(i)) + "dB   Fine";
        title(title_string);

        text(0,0.8,"Missed=" + num2str(perc_packets_missed_fine(i)));
    end
end

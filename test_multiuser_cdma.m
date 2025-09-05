%% test_multiuser_cdma.m
% Compare CDMA performance for K_users_per_slot = {1,2,4}

clear; clc; close all;

% ----- Base config -----
baseArgs = { ...
    'scheme','CDMA', ...
    'N',8, ...
    'L',128, ...
    'M',16, ...                    % (unused in CDMA path; kept for consistency)
    'bits_per_slot',128, ...
    'frames',100, ...              % bump for smoother curves
    'SNRdB_vec', -20:5:5, ...
    'threshold',0.8, ...
    'accept_ratio',0.8, ...
    'use_equal_power',true, ...    % set false + use_sic=true to demo near-far + SIC
    'use_sic',false, ...
    'compute_impostor_far',false, ...
    'compute_replay_far',false, ...
    'replay_reuse_prob',0.2 ...
};

Kset = [1 2 4];
R = cell(numel(Kset),1);

for ii = 1:numel(Kset)
    fprintf('\n=== Running K_users_per_slot = %d ===\n', Kset(ii));
    R{ii} = run_star_tdma(baseArgs{:}, 'K_users_per_slot', Kset(ii));
end

snr = R{1}.config.SNRdB_vec;

% ----- BER -----
figure; grid on;
for ii = 1:numel(Kset)
    semilogy(snr, R{ii}.BER, '-o', 'LineWidth',1.6, ...
        'DisplayName', sprintf('%d users/slot', Kset(ii)));
    hold on;
end
grid on;
xlabel('SNR (dB)'); ylabel('BER');
title('CDMA: BER vs SNR for K simultaneous users'); legend('Location','best');

% ----- FAR / FRR (legit) -----
figure; hold on; grid on;
for ii = 1:numel(Kset)
    plot(snr, R{ii}.FAR, '-s', 'LineWidth',1.6, ...
        'DisplayName', sprintf('FAR – %d users', Kset(ii)));
end
for ii = 1:numel(Kset)
    plot(snr, R{ii}.FRR, '--^', 'LineWidth',1.6, ...
        'DisplayName', sprintf('FRR – %d users', Kset(ii)));
end
xlabel('SNR (dB)'); ylabel('Rate (0..1)');
title('CDMA: FAR/FRR vs SNR'); legend('Location','best');

% ----- True impostor FAR -----
figure; hold on; grid on;
for ii = 1:numel(Kset)
    plot(snr, R{ii}.FAR_impostor, '-d', 'LineWidth',1.6, ...
        'DisplayName', sprintf('%d users', Kset(ii)));
end
xlabel('SNR (dB)'); ylabel('True Impostor FAR (0..1)');
title('CDMA: True Impostor FAR vs SNR'); legend('Location','best');

% ----- Replay FAR (conditional) -----
figure; hold on; grid on;
for ii = 1:numel(Kset)
    plot(snr, R{ii}.FAR_replay, '-o', 'LineWidth',1.6, ...
        'DisplayName', sprintf('%d users', Kset(ii)));
end
xlabel('SNR (dB)'); ylabel('Replay FAR (conditional)');
title('CDMA: Replay FAR vs SNR'); legend('Location','best');

% ----- Overall Replay Risk -----
figure; hold on; grid on;
for ii = 1:numel(Kset)
    plot(snr, R{ii}.ReplayOverall, '-o', 'LineWidth',1.6, ...
        'DisplayName', sprintf('%d users', Kset(ii)));
end
xlabel('SNR (dB)'); ylabel('Overall Replay Risk (0..1)');
title('CDMA: Overall Replay Risk vs SNR'); legend('Location','best');
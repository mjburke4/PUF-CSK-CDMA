function plot_cdma_impostors(varargin)
% plot_cdma_impostors  One-figure view of True-Impostor FAR and Replay FAR for CDMA (N=8).
%
% Usage (defaults shown):
%   plot_cdma_impostors();  % uses SNRdB_vec=-15:5:10, frames=300, reuse=0.2
%
%   % Custom sweep:
%   plot_cdma_impostors('SNRdB_vec',-25:5:10,'frames',500,'accept_ratio',0.75, ...
%       'threshold',0.7,'replay_reuse_prob',0.4,'impostor_snr_offset_db',-3, ...
%       'replay_capture_offset_db',-6);
%
% Requires: run_star_tdma.m (your current version)

% --- Defaults you can override via name-value pairs
p = inputParser;
addParameter(p,'SNRdB_vec', -15:5:10);
addParameter(p,'frames', 300);
addParameter(p,'N', 8);
addParameter(p,'L', 128);
addParameter(p,'bits_per_slot', 128);
addParameter(p,'threshold', 0.7);
addParameter(p,'accept_ratio', 0.8);
addParameter(p,'replay_reuse_prob', 0.2);      % P(reuse C) per slot
addParameter(p,'impostor_snr_db', NaN);        % fixed impostor SNR (NaN => offset from node SNR)
addParameter(p,'impostor_snr_offset_db', -3);  % impostor link SNR = node SNR + offset
addParameter(p,'replay_capture_snr_db', NaN);  % fixed capture SNR (NaN => offset from node SNR)
addParameter(p,'replay_capture_offset_db', -6);
addParameter(p,'use_node_variation', false);   % set true to randomize node SNRs
addParameter(p,'node_snr_sigma_db', 8);
parse(p, varargin{:});
cfg = p.Results;

% --- Run the star-TDMA sim (CDMA only)
res = run_star_tdma( ...
    'scheme','CDMA', ...
    'N', cfg.N, ...
    'L', cfg.L, ...
    'bits_per_slot', cfg.bits_per_slot, ...
    'frames', cfg.frames, ...
    'SNRdB_vec', cfg.SNRdB_vec, ...
    'threshold', cfg.threshold, ...
    'accept_ratio', cfg.accept_ratio, ...
    'compute_impostor_far', true, ...
    'compute_replay_far',   true, ...
    'replay_reuse_prob',    cfg.replay_reuse_prob, ...
    'impostor_snr_db',      cfg.impostor_snr_db, ...
    'impostor_snr_offset_db', cfg.impostor_snr_offset_db, ...
    'replay_capture_snr_db',  cfg.replay_capture_snr_db, ...
    'replay_capture_offset_db', cfg.replay_capture_offset_db, ...
    'use_node_variation',   cfg.use_node_variation, ...
    'node_snr_sigma_db',    cfg.node_snr_sigma_db, ...
    'K_users_per_slot', 1, ...        % single CDMA user per slot for a clean baseline
    'use_equal_power', true, ...
    'use_sic', false ...
    );

snr = res.config.SNRdB_vec;

% --- One plot: True Impostor FAR & Replay FAR (+ optional overall risk)
figure; hold on; grid on;
plot(snr, res.FAR_impostor, '-d', 'LineWidth', 1.8, 'DisplayName', 'True Impostor FAR');
plot(snr, res.FAR_replay,   '-o', 'LineWidth', 1.8, 'DisplayName', 'Replay FAR (conditional)');
plot(snr, res.ReplayOverall,'--',  'LineWidth', 1.8, 'DisplayName', 'Overall Replay Risk');

xlabel('Mean SNR (dB)');
ylabel('Rate (0..1)');
title(sprintf('CDMA Adversary FAR vs SNR  (N=%d, L=%d, bits/slot=%d, acc=%.2f, thr=%.2f, reuse=%.2f)', ...
    res.config.N, res.config.L, res.config.bits_per_slot, res.config.accept_ratio, res.config.threshold, res.config.replay_reuse_prob));
legend('Location','best');

% Also print a quick text summary in the console
fprintf('\n=== Summary ===\n');
fprintf('SNR(dB) : %s\n', mat2str(snr));
fprintf('True Impostor FAR   : %s\n', mat2str(res.FAR_impostor,3));
fprintf('Replay FAR (cond.)  : %s\n', mat2str(res.FAR_replay,3));
fprintf('Overall Replay Risk : %s\n\n', mat2str(res.ReplayOverall,3));
end


function plot_attack_success_vs_reuse_auth(varargin)
% PLOT_ATTACK_SUCCESS_VS_REUSE_AUTH
% Figure #3 (replay): Attack success vs nonce reuse, using run_star_tdma_auth.m
%
% This plots:
%   • FAR_replay (conditional acceptance rate when a replay is attempted)
%   • ReplayOverall (overall compromise rate per slot)
% Optionally overlays naive impostor FAR.
%
% Usage:
%   plot_attack_success_vs_reuse_auth('reuse_probs',0:0.1:0.9, 'SNRdB',8, ...
%       'K',2,'L',128,'codebook','walsh-orth', 'frames',30, ...
%       'threshold',0.80,'accept_ratio',0.80,'use_sic',false, 'showImpostor',true);
%
% Dependencies: run_star_tdma_auth.m on path.

p = inputParser;
addParameter(p,'reuse_probs', 0:0.1:0.9);
addParameter(p,'SNRdB', 8);
addParameter(p,'K', 2);
addParameter(p,'L', 128);
addParameter(p,'codebook', 'walsh-orth');
addParameter(p,'frames', 30);
addParameter(p,'bits_per_slot', 128);
addParameter(p,'use_sic', false);
addParameter(p,'use_equal_power', true);
addParameter(p,'rng_seed', 1337);
addParameter(p,'threshold', 0.80);
addParameter(p,'accept_ratio', 0.80);
addParameter(p,'showImpostor', true);
parse(p, varargin{:});
cfg = p.Results;

reuse_probs = cfg.reuse_probs(:).';
nr = numel(reuse_probs);

FARreplay   = nan(1, nr);
OverallRisk = nan(1, nr);
FARimp      = nan(1, nr);

fprintf('Running replay success vs nonce reuse (AUTH variant)...\n');
for ir = 1:nr
    p_reuse = reuse_probs(ir);
    res = run_star_tdma_auth( ...
        'scheme','CDMA', ...
        'N', max(8, 2*cfg.K), ...
        'K_users_per_slot', cfg.K, ...
        'cdma_codebook', cfg.codebook, ...
        'use_slot_dither', true, ...
        'use_equal_power', cfg.use_equal_power, ...
        'use_sic', cfg.use_sic, ...
        'L', cfg.L, ...
        'bits_per_slot', cfg.bits_per_slot, ...
        'frames', cfg.frames, ...
        'SNRdB_vec', cfg.SNRdB, ...
        'threshold', cfg.threshold, ...
        'accept_ratio', cfg.accept_ratio, ...
        'use_nonce', true, ...
        'replay_reuse_prob', p_reuse, ...
        'compute_replay_far', true, ...
        'compute_impostor_far', cfg.showImpostor, ...
        'log_level', 0, ...
        'rng_seed', cfg.rng_seed);

    FARreplay(ir)   = res.FAR_replay(1);
    OverallRisk(ir) = res.ReplayOverall(1);
    if cfg.showImpostor, FARimp(ir) = res.FAR_impostor(1); end
end

figure('Name','Attack success vs nonce reuse (AUTH)','Color','w'); hold on; grid on;
plot(reuse_probs, FARreplay, '-o', 'LineWidth', 1.8, 'DisplayName', 'Replay success (conditional)');
plot(reuse_probs, OverallRisk, '--s', 'LineWidth', 1.8, 'DisplayName', 'Overall replay risk');
if cfg.showImpostor
    plot(reuse_probs, FARimp, ':^', 'LineWidth', 1.2, 'DisplayName', 'Naive impostor FAR');
end
xlabel('Nonce reuse probability p'); ylabel('Attack success rate');
title(sprintf('Replay success vs reuse (SNR=%g dB, K=%d, N=%d, %s, \\tau=%.2f, acc=%.2f)', ...
    cfg.SNRdB, cfg.K, cfg.L, cfg.codebook, cfg.threshold, cfg.accept_ratio), 'Interpreter','none');
legend('Location','northwest'); legend boxoff;
ylim([0 1]);
end

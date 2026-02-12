
function plot_replay_quad_vs_reuse_auth_pairs(varargin)
% PLOT_REPLAY_QUAD_VS_REUSE_AUTH_PAIRS
% Quad plot (2x2) of replay success vs nonce reuse using run_star_tdma_auth.m,
% with FOUR explicitly specified (threshold, accept_ratio) pairs:
%   Panel 1 -> (ThrVals(1), AccVals(1))
%   Panel 2 -> (ThrVals(2), AccVals(2))
%   Panel 3 -> (ThrVals(3), AccVals(3))
%   Panel 4 -> (ThrVals(4), AccVals(4))
%
% Usage (example):
%   plot_replay_quad_vs_reuse_auth_pairs('reuse_probs', 0:0.1:0.9, ...
%       'SNRdB', 8, 'K', 2, 'L', 128, 'codebook', 'walsh-orth', 'frames', 30, ...
%       'bits_per_slot', 128, 'use_sic', false, 'use_equal_power', true, ...
%       'ThrVals', [0.70 0.75 0.80 0.85], 'AccVals', [0.70 0.80 0.85 0.90], ...
%       'showImpostor', true, 'rng_seed', 1337);
%
% Name-Value options:
%   'reuse_probs'     : vector of reuse probabilities p (default 0:0.1:0.9)
%   'SNRdB'           : SNR (dB) at BS (default 8)
%   'K'               : concurrent users per slot (default 2)
%   'L'               : spreading length (default 128)
%   'codebook'        : 'walsh-orth' | 'puf-random' (default 'walsh-orth')
%   'frames'          : Monte Carlo frames per SNR point (default 30)
%   'bits_per_slot'   : payload bits per slot (default 128)
%   'use_sic'         : enable SIC (default false)
%   'use_equal_power' : equal-power users (default true)
%   'rng_seed'        : RNG seed (default 1337)
%   'ThrVals'         : 1x4 vector of thresholds τ (default [0.75 0.80 0.85 0.90])
%   'AccVals'         : 1x4 vector of accept_ratio values (default [0.75 0.80 0.85 0.90])
%   'showImpostor'    : overlay naive impostor FAR (default true)
%
% Outputs: A 2x2 tiled figure with replay success (conditional) and overall replay risk.
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
addParameter(p,'ThrVals', [0.75 0.80 0.85 0.90]);
addParameter(p,'AccVals', [0.75 0.80 0.85 0.90]);
addParameter(p,'showImpostor', true);
parse(p, varargin{:});
cfg = p.Results;

assert(numel(cfg.ThrVals)==4 && numel(cfg.AccVals)==4, 'ThrVals and AccVals must each have 4 entries.');

reuse_probs = cfg.reuse_probs(:).';

figure('Name','Replay vs reuse (quad; explicit pairs, AUTH)','Color','w');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% Colors/markers
c1 = [0.00 0.45 0.74]; % blue
c2 = [0.85 0.33 0.10]; % orange
mk = {'o','s','^','d'};

for panel = 1:4
    thr = cfg.ThrVals(panel);
    acc = cfg.AccVals(panel);

    FARreplay   = nan(size(reuse_probs));
    OverallRisk = nan(size(reuse_probs));
    FARimp      = nan(size(reuse_probs));

    for ir = 1:numel(reuse_probs)
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
            'threshold', thr, ...
            'accept_ratio', acc, ...
            'use_nonce', true, ...
            'replay_reuse_prob', p_reuse, ...
            'compute_replay_far', true, ...
            'compute_impostor_far', cfg.showImpostor, ...
            'log_level', 0, ...
            'rng_seed', cfg.rng_seed + panel);  % small variation per panel to de-correlate

        FARreplay(ir)   = res.FAR_replay(1);
        OverallRisk(ir) = res.ReplayOverall(1);
        if cfg.showImpostor, FARimp(ir) = res.FAR_impostor(1); end
    end

    nexttile; hold on; grid on; ylim([0 1]);
    plot(reuse_probs, FARreplay, ['-' mk{1}], 'Color', c1, 'LineWidth', 1.8, ...
         'DisplayName','Replay success (cond.)');
    plot(reuse_probs, OverallRisk, ['--' mk{2}], 'Color', c2, 'LineWidth', 1.8, ...
         'DisplayName','Overall replay risk');
    if cfg.showImpostor
        plot(reuse_probs, FARimp, ':', 'Color', 0.5*[1 1 1], 'LineWidth', 1.2, ...
             'DisplayName','Naive impostor FAR');
    end

    xlabel('Nonce reuse probability p'); ylabel('Attack success rate');
    title(sprintf('\\tau=%.2f, acc=%.2f  |  SNR=%g dB, K=%d, N=%d, %s', ...
          thr, acc, cfg.SNRdB, cfg.K, cfg.L, cfg.codebook), 'Interpreter','none');

    if panel==1
        legend('Location','northwest'); legend boxoff;
    end
end
end

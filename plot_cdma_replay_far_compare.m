
function plot_cdma_replay_far_compare(varargin)
% plot_cdma_replay_far_compare
% Compare replay attacker effectiveness for three CDMA code strategies:
%   1) PUF-random spreading            (unpredictable)
%   2) Walsh + PUF mask (slot dither)  (orthogonal + unpredictable across slots)
%   3) Walsh-only                      (orthogonal, predictable)
%
% Plots two curves per mode:
%   - Conditional Replay FAR (given a replay trial occurred)
%   - Overall Replay Risk    (accepted replays / total slots)
%
% IMPORTANT: With proper freshness (nonce + no reuse), there will be ~no trials.
% To visualize conditional FAR, we induce rare pair reuse via 'replay_reuse_prob' (e.g., 0.01).
%
% Optional Name-Value args:
%   'K'                   (default 4)       : simultaneous users per slot
%   'L'                   (default 128)     : spreading length (power of two)
%   'frames'              (default 30)      : Monte Carlo frames per SNR
%   'bits_per_slot'       (default 128)
%   'SNRdB_vec'           (default -10:2:12)
%   'use_sic'             (default false)
%   'use_equal_power'     (default true)
%   'rng_seed'            (default 1337)
%   'replay_reuse_prob'   (default 0.01)    : small reuse to generate replay trials
%   'replay_capture_off'  (default -6)      : capture SNR offset from node SNR (dB)
%   'impostor_snr_off'    (default -3)      : not used here; keep for completeness
%   'use_channel_factor'  (default false)   : enable channel similarity gate
%   'tau_h'               (default 0.75)
%   'coherence_prob'      (default 1.0)     : prob the channel CHANGES between capture and replay
%
% Example:
%   plot_cdma_replay_far_compare('K',4,'use_channel_factor',true,'tau_h',0.8);

    % ---------------- Parse inputs
    p = inputParser;
    addParameter(p,'K',4);
    addParameter(p,'L',128);
    addParameter(p,'frames',30);
    addParameter(p,'bits_per_slot',128);
    addParameter(p,'SNRdB_vec',-25:5:5);
    addParameter(p,'use_sic',false);
    addParameter(p,'use_equal_power',true);
    addParameter(p,'rng_seed',1337);
    addParameter(p,'replay_reuse_prob',0.01);
    addParameter(p,'replay_capture_off',-6);
    addParameter(p,'impostor_snr_off',-3);
    addParameter(p,'use_channel_factor',false);
    addParameter(p,'tau_h',0.75);
    addParameter(p,'coherence_prob',1.0);
    parse(p,varargin{:});
    cfg = p.Results;

    assert(abs(log2(cfg.L) - round(log2(cfg.L))) < 1e-9, 'L must be a power of two.');

    % ---------------- Modes
    modes = { ...
        struct('name','PUF-random',      'cdma_codebook','puf-random', 'use_slot_dither',true),  ...
        struct('name','Walsh+PUF-mask',  'cdma_codebook','walsh-orth', 'use_slot_dither',true), ...
        struct('name','Walsh-only',      'cdma_codebook','walsh-orth', 'use_slot_dither',false) ...
    };
    nModes = numel(modes);
    numSNR = numel(cfg.SNRdB_vec);

    FARc = nan(nModes, numSNR);   % conditional replay FAR
    Rrisk= nan(nModes, numSNR);   % overall replay risk
    Rtr  = nan(nModes, numSNR);   % number of replay trials (for sanity)

    % ---------------- Run sweeps
    fprintf('Running CDMA replay FAR sweeps (K=%d, reuse_prob=%.3f)...\n', cfg.K, cfg.replay_reuse_prob);
    for m = 1:nModes
        M = modes{m};
        fprintf('  Mode: %s\n', M.name);
        try
            res = run_star_tdma( ...
                'scheme','CDMA', ...
                'N', max(6, 2*cfg.K), ...
                'K_users_per_slot', cfg.K, ...
                'cdma_codebook', M.cdma_codebook, ...
                'use_slot_dither', M.use_slot_dither, ...
                'use_equal_power', cfg.use_equal_power, ...
                'use_sic', cfg.use_sic, ...
                'L', cfg.L, ...
                'bits_per_slot', cfg.bits_per_slot, ...
                'frames', cfg.frames, ...
                'SNRdB_vec', cfg.SNRdB_vec, ...
                'use_nonce', true, ...
                'replay_reuse_prob', cfg.replay_reuse_prob, ...
                'replay_capture_snr_db', NaN, ...
                'replay_capture_offset_db', cfg.replay_capture_off, ...
                'impostor_snr_db', NaN, ...
                'impostor_snr_offset_db', cfg.impostor_snr_off, ...
                'compute_replay_far', true, ...
                'compute_impostor_far', false, ...   % not needed here
                'use_channel_factor', cfg.use_channel_factor, ...
                'tau_h', cfg.tau_h, ...
                'coherence_prob', cfg.coherence_prob, ...
                'log_level', 0, ...
                'rng_seed', cfg.rng_seed);
            FARc(m, :) = res.FAR_replay;
            Rrisk(m, :)= res.ReplayOverall;
            Rtr(m, :)  = res.replay_trials;
        catch ME
            warning('Mode %s failed: %s', M.name, ME.message);
        end
    end

    % ---------------- Plot Conditional FAR
    figure; hold on; grid on;
    title(sprintf('Conditional Replay FAR vs SNR (K=%d, L=%d, frames=%d, reuse=%.3f, CF=%d)', ...
        cfg.K, cfg.L, cfg.frames, cfg.replay_reuse_prob, cfg.use_channel_factor));
    xlabel('SNR (dB)'); ylabel('Replay FAR (conditional)');
    styles = {'-o','--s',':^'};
    for m = 1:nModes
        plot(cfg.SNRdB_vec, FARc(m,:), styles{m}, 'LineWidth', 1.7, 'DisplayName', modes{m}.name);
    end
    ylim([0 1]); legend('Location','northwest');

    % ---------------- Plot Overall Replay Risk
    figure; hold on; grid on;
    title(sprintf('Overall Replay Risk vs SNR (K=%d, L=%d, frames=%d, reuse=%.3f, CF=%d)', ...
        cfg.K, cfg.L, cfg.frames, cfg.replay_reuse_prob, cfg.use_channel_factor));
    xlabel('SNR (dB)'); ylabel('Overall Replay Compromise Rate');
    for m = 1:nModes
        plot(cfg.SNRdB_vec, Rrisk(m,:), styles{m}, 'LineWidth', 1.7, 'DisplayName', modes{m}.name);
    end
    ylim([0 1]); legend('Location','northwest');

    % ---------------- Optional console sanity: how many trials?
    fprintf('\nReplay trials per SNR (rows=modes, cols=SNR):\n');
    disp(Rtr);
end

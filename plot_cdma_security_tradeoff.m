
function plot_cdma_security_tradeoff(varargin)
% plot_cdma_security_tradeoff
% Compare CDMA security/performance across three code strategies:
%   1) PUF-random spreading            -> unpredictable, non-orthogonal
%   2) Walsh + PUF mask (slot dither)  -> orthogonal within-slot + unpredictable across slots
%   3) Walsh-only (no mask)            -> orthogonal but predictable (no secrecy)
%
% Outputs two figures:
%   - BER vs SNR (log scale)
%   - FRR vs SNR (linear)
%
% Requirements: run_star_tdma.m on path (version with 'cdma_codebook' and 'use_slot_dither').
%
% Optional Name-Value args:
%   'K'               (default 4)     : simultaneous CDMA users per slot
%   'L'               (default 128)   : spreading length (power of two)
%   'frames'          (default 24)    : Monte Carlo frames per SNR
%   'bits_per_slot'   (default 128)   : payload bits per slot
%   'SNRdB_vec'       (default -10:2:12)
%   'use_sic'         (default false)
%   'use_equal_power' (default true)
%   'rng_seed'        (default 1337)
%
% Example:
%   plot_cdma_security_tradeoff('K',4,'L',128,'frames',30,'SNRdB_vec',-8:2:12,'use_sic',true);

    % ---------------- Parse inputs
    p = inputParser;
    addParameter(p,'K',4);
    addParameter(p,'L',128);
    addParameter(p,'frames',24);
    addParameter(p,'bits_per_slot',128);
    addParameter(p,'SNRdB_vec',-10:2:12);
    addParameter(p,'use_sic',false);
    addParameter(p,'use_equal_power',true);
    addParameter(p,'rng_seed',1337);
    parse(p,varargin{:});
    cfg = p.Results;

    assert(abs(log2(cfg.L) - round(log2(cfg.L))) < 1e-9, 'L must be a power of two.');

    % ---------------- Configurations to compare
    modes = { ...
        struct('name','PUF-random',        'cdma_codebook','puf-random', 'use_slot_dither',true),  ...
        struct('name','Walsh+PUF-mask',    'cdma_codebook','walsh-orth', 'use_slot_dither',true), ...
        struct('name','Walsh-only',        'cdma_codebook','walsh-orth', 'use_slot_dither',false) ...
    };
    nModes = numel(modes);
    numSNR = numel(cfg.SNRdB_vec);

    BER = nan(nModes, numSNR);
    FRR = nan(nModes, numSNR);

    % ---------------- Run sweeps
    fprintf('Running CDMA sweeps for K=%d users...\n', cfg.K);
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
                'replay_reuse_prob', 0.0, ...      % freshness on; replay irrelevant here
                'compute_replay_far', false, ...
                'compute_impostor_far', false, ...
                'log_level', 0, ...
                'rng_seed', cfg.rng_seed);
            BER(m, :) = res.BER;
            FRR(m, :) = res.FRR;
        catch ME
            warning('Mode %s failed: %s', M.name, ME.message);
        end
    end

    % ---------------- Plot BER
    figure; hold on; grid on;
    title(sprintf('CDMA BER vs SNR (K=%d, L=%d, frames=%d, bits/slot=%d, SIC=%d)', ...
        cfg.K, cfg.L, cfg.frames, cfg.bits_per_slot, cfg.use_sic));
    xlabel('SNR (dB)'); ylabel('BER'); set(gca,'YScale','log');
    styles = {'-o','--s',':^'}; % one style per mode
    for m = 1:nModes
        plot(cfg.SNRdB_vec, BER(m,:), styles{m}, 'LineWidth', 1.7, 'DisplayName', modes{m}.name);
    end
    legend('Location','southwest');

    % ---------------- Plot FRR
    figure; hold on; grid on;
    title(sprintf('CDMA FRR vs SNR (K=%d, L=%d, frames=%d, bits/slot=%d, SIC=%d)', ...
        cfg.K, cfg.L, cfg.frames, cfg.bits_per_slot, cfg.use_sic));
    xlabel('SNR (dB)'); ylabel('FRR (0..1)');
    for m = 1:nModes
        plot(cfg.SNRdB_vec, FRR(m,:), styles{m}, 'LineWidth', 1.7, 'DisplayName', modes{m}.name);
    end
    ylim([0 1]); legend('Location','northeast');
end

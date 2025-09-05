
function plot_cdma_sic_frr_effect(varargin)
% plot_cdma_sic_frr_effect
% Evaluate the impact of Successive Interference Cancellation (SIC) on FRR
% in your CDMA simulator (run_star_tdma.m).
%
% It compares FRR vs SNR for:
%   - K ∈ {2, 4} simultaneous users (default; configurable)
%   - Power cases: default equal-power (sigma=0). You can add near-far sigma, e.g., [0 6 12].
% and overlays results for SIC OFF vs SIC ON.
%
% It also prints a table of average FRR across the SNR sweep to quantify differences.
%
% Optional Name-Value args:
%   'L'               (default 128)
%   'frames'          (default 24)
%   'bits_per_slot'   (default 128)
%   'SNRdB_vec'       (default -8:2:12)
%   'K_list'          (default [2 4])
%   'sigma_list'      (default 0)      % 0 => equal power; add 6 or 12 for near-far
%   'accept_ratio'    (default 0.8)    % keep consistent with your runs
%   'threshold'       (default 0.8)
%   'rng_seed'        (default 1337)
%
% Example:
%   plot_cdma_sic_frr_effect('L',128,'frames',30,'K_list',[2 4],'sigma_list',0);
%
% Notes:
% - Uses 'cdma_codebook' = 'puf-random' to isolate SIC effects (no Walsh changes).
% - Replay/impostor are disabled for speed.
% - Nonce is ON (freshness), but irrelevant to FRR (legit decision path).

    % ---------------- Parse inputs
    p = inputParser;
    addParameter(p,'L',128);
    addParameter(p,'frames',24);
    addParameter(p,'bits_per_slot',128);
    addParameter(p,'SNRdB_vec',-8:2:12);
    addParameter(p,'K_list',[2 4]);
    addParameter(p,'sigma_list',0);
    addParameter(p,'accept_ratio',0.8);
    addParameter(p,'threshold',0.8);
    addParameter(p,'rng_seed',1337);
    parse(p,varargin{:});
    cfg = p.Results;

    assert(abs(log2(cfg.L) - round(log2(cfg.L))) < 1e-9, 'L must be a power of two.');

    Ks = cfg.K_list(:)';
    sigmas = cfg.sigma_list(:)';
    if isempty(sigmas), sigmas = 0; end

    % Storage: results{Kindex, sigmaIndex, sicFlag}
    results = cell(numel(Ks), numel(sigmas), 2); % 1: SIC off, 2: SIC on

    fprintf('Running FRR SIC-effect sweeps...\n');
    for iK = 1:numel(Ks)
        K = Ks(iK);
        for isg = 1:numel(sigmas)
            sigma_db = sigmas(isg);
            for sicFlag = 0:1
                fprintf('  K=%d, sigma=%2d dB, SIC=%d ... ', K, sigma_db, sicFlag);
                try
                    res = run_star_tdma( ...
                        'scheme','CDMA', ...
                        'N', max(6, 2*K), ...
                        'K_users_per_slot', K, ...
                        'cdma_codebook', 'puf-random', ...
                        'use_slot_dither', true, ...
                        'use_equal_power', (sigma_db==0), ...
                        'power_sigma_db', max(1e-9, sigma_db), ...
                        'use_sic', logical(sicFlag), ...
                        'L', cfg.L, ...
                        'bits_per_slot', cfg.bits_per_slot, ...
                        'frames', cfg.frames, ...
                        'SNRdB_vec', cfg.SNRdB_vec, ...
                        'accept_ratio', cfg.accept_ratio, ...
                        'threshold', cfg.threshold, ...
                        'use_nonce', true, ...
                        'replay_reuse_prob', 0.0, ...
                        'compute_replay_far', false, ...
                        'compute_impostor_far', false, ...
                        'log_level', 0, ...
                        'rng_seed', cfg.rng_seed);
                    results{iK, isg, sicFlag+1} = res;
                    fprintf('done.\n');
                catch ME
                    fprintf('ERROR: %s\n', ME.message);
                    results{iK, isg, sicFlag+1} = [];
                end
            end
        end
    end

    % Plot FRR overlays per K and sigma
    colors = lines(3);
    for iK = 1:numel(Ks)
        K = Ks(iK);
        figure; tiledlayout(1, numel(sigmas),'Padding','compact','TileSpacing','compact');
        for isg = 1:numel(sigmas)
            sigma_db = sigmas(isg);
            nexttile; hold on; grid on;
            title(sprintf('K=%d, %s', K, sigma_label(sigma_db)));
            xlabel('SNR (dB)'); ylabel('FRR (0..1)');
            % SIC OFF
            r0 = results{iK,isg,1};
            if ~isempty(r0)
                plot(cfg.SNRdB_vec, r0.FRR, '-o', 'LineWidth', 1.6, 'DisplayName', 'SIC off', 'Color', colors(1,:));
            end
            % SIC ON
            r1 = results{iK,isg,2};
            if isempty(r1) % <- MATLAB '!' not allowed, fix immediately below
                plot(cfg.SNRdB_vec, r1.FRR, '--s', 'LineWidth', 1.6, 'DisplayName', 'SIC on', 'Color', colors(2,:));
            end
            % correct MATLAB syntax:
            if ~isempty(r1)
                plot(cfg.SNRdB_vec, r1.FRR, '--s', 'LineWidth', 1.6, 'DisplayName', 'SIC on', 'Color', colors(2,:));
            end
            ylim([0 1]); legend('Location','southwest');
        end
        sgtitle(sprintf('CDMA FRR vs SNR — SIC effect (L=%d, frames=%d, bits/slot=%d, thr=%.2f, acc=%.2f)', ...
            cfg.L, cfg.frames, cfg.bits_per_slot, cfg.threshold, cfg.accept_ratio));
    end

    % Print average FRR table across SNRs
    fprintf('\nAverage FRR across SNR sweep:\n');
    fprintf('K  sigma(dB)   SIC   avgFRR\n');
    fprintf('-- ---------   ---   ------\n');
    for iK = 1:numel(Ks)
        K = Ks(iK);
        for isg = 1:numel(sigmas)
            sigma_db = sigmas(isg);
            for sicFlag = 0:1
                r = results{iK,isg,sicFlag+1};
                if isempty(r), avgFRR = NaN; else, avgFRR = mean(r.FRR); end
                fprintf('%d     %2d       %d    %8.4g\n', K, sigma_db, sicFlag, avgFRR);
            end
        end
    end

    if nargout > 0
        varargout{1} = results;
    end
end

function s = sigma_label(sigma_db)
    if sigma_db==0, s = 'equal power';
    else, s = sprintf('near-far (\\sigma=%d dB)', sigma_db);
    end
end

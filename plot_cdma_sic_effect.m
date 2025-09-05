
function plot_cdma_sic_effect(varargin)
% plot_cdma_sic_effect
% Evaluate the impact of Successive Interference Cancellation (SIC) on CDMA
% in your run_star_tdma simulator.
%
% It compares BER vs SNR for:
%   - K ∈ {2, 4} simultaneous users
%   - Power distributions: equal-power and near-far with sigma ∈ {6, 12} dB
% and overlays results for SIC OFF vs SIC ON.
%
% It also prints a small table of average BER across the SNR sweep to quantify gains.
%
% Requirements: run_star_tdma.m on path (version with multi-user CDMA + SIC toggle).
%
% Optional Name-Value args:
%   'L'               (default 128)
%   'frames'          (default 24)
%   'bits_per_slot'   (default 128)
%   'SNRdB_vec'       (default -8:2:12)
%   'rng_seed'        (default 1337)
%
% Example:
%   plot_cdma_sic_effect('L',128,'frames',30,'SNRdB_vec',-8:2:12);

    % ---------------- Parse inputs
    p = inputParser;
    addParameter(p,'L',128);
    addParameter(p,'frames',24);
    addParameter(p,'bits_per_slot',128);
    addParameter(p,'SNRdB_vec',-8:2:12);
    addParameter(p,'rng_seed',1337);
    parse(p,varargin{:});
    cfg = p.Results;

    assert(abs(log2(cfg.L) - round(log2(cfg.L))) < 1e-9, 'L must be a power of two.');

    % ---------------- Experiment matrix
    Ks = [2 4];
    pow_sigmas = [0 6 12];  % 0 => equal power; otherwise near-far sigma (dB)

    % Preallocate storage: results{Kindex, sigmaIndex, sicFlag} -> struct
    results = cell(numel(Ks), numel(pow_sigmas), 2); % 1: SIC off, 2: SIC on

    fprintf('Running SIC effect sweeps...\n');
    for iK = 1:numel(Ks)
        K = Ks(iK);
        for isg = 1:numel(pow_sigmas)
            sigma_db = pow_sigmas(isg);
            for sicFlag = 0:1
                fprintf('  K=%d, sigma=%2d dB, SIC=%d ... ', K, sigma_db, sicFlag);
                try
                    res = run_star_tdma( ...
                        'scheme','CDMA', ...
                        'N', max(6, 2*K), ...
                        'K_users_per_slot', K, ...
                        'cdma_codebook', 'puf-random', ...   % focus on SIC, not Walsh
                        'use_slot_dither', true, ...
                        'use_equal_power', (sigma_db==0), ...
                        'power_sigma_db', max(1e-9, sigma_db), ...
                        'use_sic', logical(sicFlag), ...
                        'L', cfg.L, ...
                        'bits_per_slot', cfg.bits_per_slot, ...
                        'frames', cfg.frames, ...
                        'SNRdB_vec', cfg.SNRdB_vec, ...
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

    % ---------------- Plot BER overlays per K and sigma
    colors = lines(3);  % just for consistency across figures
    for iK = 1:numel(Ks)
        K = Ks(iK);
        figure; tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
        for isg = 1:numel(pow_sigmas)
            sigma_db = pow_sigmas(isg);
            nexttile; hold on; grid on;
            title(sprintf('K=%d, %s', K, sigma_label(sigma_db)));
            xlabel('SNR (dB)'); ylabel('BER'); set(gca,'YScale','log');
            % SIC OFF
            r0 = results{iK,isg,1};
            if ~isempty(r0)
                plot(cfg.SNRdB_vec, r0.BER, '-o', 'LineWidth', 1.6, 'DisplayName', 'SIC off', 'Color', colors(1,:));
            end
            % SIC ON
            r1 = results{iK,isg,2};
            if ~isempty(r1)
                plot(cfg.SNRdB_vec, r1.BER, '--s', 'LineWidth', 1.6, 'DisplayName', 'SIC on', 'Color', colors(2,:));
            end
            legend('Location','southwest');
        end
        sgtitle(sprintf('CDMA BER vs SNR — SIC effect (L=%d, frames=%d, bits/slot=%d)', cfg.L, cfg.frames, cfg.bits_per_slot));
    end

    % ---------------- Print average BER table across SNRs
    fprintf('\nAverage BER across SNR sweep:\n');
    fprintf('K  sigma(dB)   SIC   avgBER\n');
    fprintf('-- ---------   ---   ------\n');
    for iK = 1:numel(Ks)
        K = Ks(iK);
        for isg = 1:numel(pow_sigmas)
            sigma_db = pow_sigmas(isg);
            for sicFlag = 0:1
                r = results{iK,isg,sicFlag+1};
                if isempty(r), avgBER = NaN; else, avgBER = mean(r.BER); end
                fprintf('%d     %2d       %d    %8.4g\n', K, sigma_db, sicFlag, avgBER);
            end
        end
    end

    % ------------- Optional: return results if caller asks
    if nargout > 0
        varargout{1} = results;
    end
end

function s = sigma_label(sigma_db)
    if sigma_db==0
        s = 'equal power';
    else
        s = sprintf('near-far (\\sigma=%d dB)', sigma_db);
    end
end

function plot_cdma_ber_compare(varargin)
% plot_cdma_ber_compare
% Compare CDMA BER vs SNR for K = [1, 2, 4] under:
%   cdma_codebook = 'puf-random'  and  'walsh-orth'
%
% Usage examples:
%   plot_cdma_ber_compare;  % defaults
%   plot_cdma_ber_compare('L',128,'frames',30,'SNRdB_vec',-8:2:12,'use_sic',true);

    % ---------------- Parse inputs
    p = inputParser;
    addParameter(p,'L',128);
    addParameter(p,'frames',20);
    addParameter(p,'bits_per_slot',128);
    addParameter(p,'SNRdB_vec',-10:2:12);
    addParameter(p,'use_sic',false);
    addParameter(p,'use_equal_power',true);
    addParameter(p,'rng_seed',1337);
    parse(p,varargin{:});
    cfg = p.Results;

    assert(abs(log2(cfg.L) - round(log2(cfg.L))) < 1e-9, 'L must be a power of two.');

    % ---------------- Experiment grid
    Ks        = [1 2 4];
    codebooks = {'puf-random','walsh-orth'};
    numK   = numel(Ks);
    numSNR = numel(cfg.SNRdB_vec);

    % line styles per codebook (no hyphens in field names needed)
    styles = { {'-o','-s','-^'}, {'--o','--s','--^'} };  % {codebook}{K-index}

    % ---------------- Storage: BER(codebook_index, K_index, snr_index)
    BER = nan(numel(codebooks), numK, numSNR);
    N = 8; % added by mike
    % ---------------- Run sweeps
    fprintf('Running CDMA BER sweeps...\n');
    for c = 1:numel(codebooks)
        cb = codebooks{c};
        fprintf('  codebook = %s\n', cb);
        for ik = 1:numK
            K = Ks(ik);
            fprintf('    K_users_per_slot = %d ... ', K);
            try
                res = run_star_tdma( ...
                    'scheme','CDMA', ...
                    'N', N, ... %max(4, 2*K), ... Changed my Mike 20250831
                    'K_users_per_slot', K, ...
                    'cdma_codebook', cb, ...
                    'use_slot_dither', true, ...
                    'use_equal_power', cfg.use_equal_power, ...
                    'use_sic', cfg.use_sic, ...
                    'L', cfg.L, ...
                    'bits_per_slot', cfg.bits_per_slot, ...
                    'frames', cfg.frames, ...
                    'SNRdB_vec', cfg.SNRdB_vec, ...
                    'use_nonce', true, ...
                    'replay_reuse_prob', 0.0, ...      % freshness; replay off for BER
                    'compute_replay_far', false, ...
                    'compute_impostor_far', false, ...
                    'log_level', 0, ...
                    'rng_seed', cfg.rng_seed);

                BER(c, ik, :) = reshape(res.BER, 1, 1, []);
                fprintf('done.\n');
            catch ME
                fprintf('ERROR: %s\n', ME.message);
            end
        end
    end

    % ---------------- Plot
    figure; hold on; grid on;
    title(sprintf('CDMA BER vs SNR — PUF-random vs Walsh-orth (L=%d, frames=%d, bits/slot=%d, SIC=%d, N=%d)', ...
        cfg.L, cfg.frames, cfg.bits_per_slot, cfg.use_sic, N));
    xlabel('SNR (dB)'); ylabel('BER');

    % Plot PUF-random first
    c = 1;
    for ik = 1:numK
        plot(cfg.SNRdB_vec, squeeze(BER(c,ik,:))', styles{c}{ik}, 'LineWidth', 1.6, ...
             'DisplayName', sprintf('PUF-rand, K=%d', Ks(ik)));
    end
    % Then Walsh-orth overlay
    c = 2;
    for ik = 1:numK
        plot(cfg.SNRdB_vec, squeeze(BER(c,ik,:))', styles{c}{ik}, 'LineWidth', 1.6, ...
             'DisplayName', sprintf('Walsh-orth, K=%d', Ks(ik)));
    end

    set(gca, 'YScale', 'log');
    legend('Location','southwest');
end
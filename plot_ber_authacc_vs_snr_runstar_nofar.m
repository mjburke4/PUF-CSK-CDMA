
function plot_ber_authacc_vs_snr_runstar_nofar(varargin)
% PLOT_BER_AUTHACC_VS_SNR_RUNSTAR_NOFAR
% Plot 2 — BER / Authentication Accuracy vs. SNR (fixed threshold, no FAR overlay) using run_star_tdma.m
%
% This is identical to plot_ber_authacc_vs_snr_runstar, except it omits the FAR overlay
% and right-hand axis for a cleaner figure suitable for the paper.
%
% Usage:
%   plot_ber_authacc_vs_snr_runstar_nofar('Lset',[64 128 256], 'Kset',[1 4], ...
%       'codebook','walsh-orth', 'SNRdB_vec',-10:2:12, 'frames',20, ...
%       'threshold',0.8, 'accept_ratio',0.8, 'use_sic',false);
%
% Dependencies: run_star_tdma.m on path.

p = inputParser;
addParameter(p,'Lset',[64 128 256]);
addParameter(p,'Kset',[1 4]);
addParameter(p,'codebook','walsh-orth');
addParameter(p,'SNRdB_vec',-10:2:12);
addParameter(p,'frames',20);
addParameter(p,'bits_per_slot',128);
addParameter(p,'use_sic',false);
addParameter(p,'use_equal_power',true);
addParameter(p,'rng_seed',1337);
addParameter(p,'threshold',0.8);
addParameter(p,'accept_ratio',0.8);
parse(p,varargin{:});
cfg = p.Results;

if ischar(cfg.codebook), codebooks = {cfg.codebook};
elseif isstring(cfg.codebook), codebooks = cellstr(cfg.codebook);
else, codebooks = cfg.codebook;
end

SNR = cfg.SNRdB_vec(:).';
colors = lines(max(numel(cfg.Lset), numel(cfg.Kset)));
markerList = {'o','s','^','d','v','>'};

for cbix = 1:numel(codebooks)
    cb = codebooks{cbix};

    BER  = struct(); FRR = struct();

    for ik = 1:numel(cfg.Kset)
        K = cfg.Kset(ik);
        kfield = sprintf('K%d',K);
        BER.(kfield) = struct(); FRR.(kfield) = struct();

        for iL = 1:numel(cfg.Lset)
            L = cfg.Lset(iL);
            fprintf('[%s] Running K=%d, L=%d ...\n', cb, K, L);
            res = run_star_tdma( ...
                'scheme','CDMA', ...
                'N', max(8, 2*K), ...
                'K_users_per_slot', K, ...
                'cdma_codebook', cb, ...
                'use_slot_dither', true, ...
                'use_equal_power', cfg.use_equal_power, ...
                'use_sic', cfg.use_sic, ...
                'L', L, ...
                'bits_per_slot', cfg.bits_per_slot, ...
                'frames', cfg.frames, ...
                'SNRdB_vec', cfg.SNRdB_vec, ...
                'threshold', cfg.threshold, ...
                'accept_ratio', cfg.accept_ratio, ...
                'use_nonce', true, ...
                'replay_reuse_prob', 0.0, ...
                'compute_replay_far', false, ...
                'compute_impostor_far', false, ...
                'log_level', 0, ...
                'rng_seed', cfg.rng_seed);

            BER.(kfield).(sprintf('L%d',L)) = res.BER(:).';
            FRR.(kfield).(sprintf('L%d',L)) = res.FRR(:).';
        end
    end

    % Figure 1: BER vs SNR
    figure('Name',sprintf('BER vs SNR (%s, fixed \\tau)', cb),'Color','w');
    tiledlayout(numel(cfg.Kset),1,'TileSpacing','compact','Padding','compact');
    for ik = 1:numel(cfg.Kset)
        K = cfg.Kset(ik); kfield = sprintf('K%d',K);
        nexttile; hold on; grid on; set(gca,'YScale','log'); ylim([1e-4 1]);
        for iL = 1:numel(cfg.Lset)
            L = cfg.Lset(iL);
            c = colors(iL,:); mk = markerList{1+mod(iL-1,numel(markerList))};
            plot(SNR, BER.(kfield).(sprintf('L%d',L)), ['-' mk], 'Color', c, 'LineWidth', 1.6, ...
                'DisplayName', sprintf('N=%d', L));
        end
        title(sprintf('BER vs SNR — K=%d (%s)', K, cb), 'Interpreter','none');
        xlabel('SNR (dB)'); ylabel('BER');
        legend('Location','southwest'); legend boxoff;
    end

    % Figure 2: Auth accuracy vs SNR (no FAR overlay)
    figure('Name',sprintf('AuthAcc vs SNR (%s, fixed \\tau)', cb),'Color','w');
    tiledlayout(numel(cfg.Kset),1,'TileSpacing','compact','Padding','compact');
    for ik = 1:numel(cfg.Kset)
        K = cfg.Kset(ik); kfield = sprintf('K%d',K);
        nexttile; hold on; grid on; ylim([0 1]);
        for iL = 1:numel(cfg.Lset)
            L = cfg.Lset(iL);
            c = colors(iL,:); mk = markerList{1+mod(iL-1,numel(markerList))};
            authacc = 1 - FRR.(kfield).(sprintf('L%d',L));
            plot(SNR, authacc, ['-' mk], 'Color', c, 'LineWidth', 1.6, ...
                 'DisplayName', sprintf('N=%d', L));
        end
        xlabel('SNR (dB)'); ylabel('Authentication accuracy (1 - FRR)');
        title(sprintf('Auth accuracy vs SNR — K=%d (%s), \\tau=%.2f, acc=%.2f', ...
              K, cb, cfg.threshold, cfg.accept_ratio), 'Interpreter','none');
        legend('Location','southeast'); legend boxoff;
    end
end
end

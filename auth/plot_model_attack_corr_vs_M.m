
function plot_model_attack_corr_vs_M(varargin)
% PLOT_MODEL_ATTACK_CORR_VS_M
% Modeling attack analysis: for each CRP leakage M, trains per-chip models,
% evaluates per-chip accuracy, slot-level attack success, and **mean correlation**
% between predicted spread and true spread (attacker closeness).
%
% Usage:
% plot_model_attack_corr_vs_M('Mset',[64 128 256 512], 'L',128, 'challenge_len',64, ...
%   'train_frames',2000, 'test_frames',400, 'SNRdB',8, 'impostor_snr_db',5, ...
%   'threshold',0.80, 'accept_ratio',0.80, 'bits_per_slot',128, ...
%   'crp_tap','expanded', 'use_seed_expansion',false, 'feature','apuf', 'rng_seed',123, ...
%   'save_csv','results_M_corr.csv');
%
% Outputs a figure with three curves vs M:
%  - attack success (slot-level)
%  - mean per-chip accuracy
%  - mean correlation (dot-product normalized) between predicted and true spread
% and optionally writes a CSV with numeric results.
%
% Dependencies: run_star_tdma_auth.m on path.

p = inputParser;
addParameter(p,'Mset', [64 128 256 512 1024]);
addParameter(p,'L', 128);
addParameter(p,'challenge_len', 64);
addParameter(p,'train_frames', 2000);
addParameter(p,'test_frames', 400);
addParameter(p,'SNRdB', 8);
addParameter(p,'impostor_snr_db', 5);
addParameter(p,'threshold', 0.80);
addParameter(p,'accept_ratio', 0.80);
addParameter(p,'bits_per_slot', 128);
addParameter(p,'crp_tap', 'expanded');           % 'expanded' | 'raw' (expanded required here)
addParameter(p,'use_seed_expansion', false);     % IMPORTANT for learnability
addParameter(p,'feature', 'apuf');               % 'apuf' | 'bits'
addParameter(p,'rng_seed', 123);
addParameter(p,'shuffle', true);
addParameter(p,'save_csv', '');
parse(p, varargin{:});
cfg = p.Results;

rng(cfg.rng_seed);

% Harvest CRPs once (train and test pools)
res_tr = run_star_tdma_auth( ...
    'scheme','CDMA', 'N', 1, 'K_users_per_slot', 1, ...
    'L', cfg.L, 'bits_per_slot', 1, ...
    'frames', cfg.train_frames, 'SNRdB_vec', 0, ...
    'use_nonce', true, 'nonce_len', cfg.challenge_len, ...
    'replay_reuse_prob', 0.0, ...
    'auth_mode','exposedCRP', 'crp_tap', cfg.crp_tap, ...
    'use_seed_expansion', cfg.use_seed_expansion, ...
    'log_level', 0, 'rng_seed', cfg.rng_seed);

res_te = run_star_tdma_auth( ...
    'scheme','CDMA', 'N', 1, 'K_users_per_slot', 1, ...
    'L', cfg.L, 'bits_per_slot', 1, ...
    'frames', cfg.test_frames, 'SNRdB_vec', 0, ...
    'use_nonce', true, 'nonce_len', cfg.challenge_len, ...
    'replay_reuse_prob', 0.0, ...
    'auth_mode','exposedCRP', 'crp_tap', cfg.crp_tap, ...
    'use_seed_expansion', cfg.use_seed_expansion, ...
    'log_level', 0, 'rng_seed', cfg.rng_seed + 1);

% Pack CRPs (expanded tap required)
if ~strcmpi(cfg.crp_tap,'expanded')
    error('plot_model_attack_corr_vs_M requires crp_tap=''expanded'' (full L-bit spreads).');
end
[Ctr, Rtr] = pack_crp_logs_expanded(res_tr.CRP_LOGS, cfg.challenge_len, cfg.L);
[Cte, Rte] = pack_crp_logs_expanded(res_te.CRP_LOGS, cfg.challenge_len, cfg.L);
Lchips = size(Rtr,2);

% Optionally shuffle training set
if cfg.shuffle
    idx_tr = randperm(size(Ctr,1));
    Ctr = Ctr(idx_tr,:);
    Rtr = Rtr(idx_tr,:);
end

% Build features
if strcmpi(cfg.feature,'apuf')
    Xtr = apuf_features(Ctr);
    Xte = apuf_features(Cte);
else
    Xtr = [ones(size(Ctr,1),1), double(Ctr)];
    Xte = [ones(size(Cte,1),1), double(Cte)];
end

Mset = cfg.Mset(:).';
nM = numel(Mset);
attack_success = nan(1,nM);
perchip_acc = nan(1,nM);
mean_corr = nan(1,nM);

for im = 1:nM
    M = min(Mset(im), size(Xtr,1));
    X_M = Xtr(1:M,:);
    Y_M = Rtr(1:M,:);  % M x L

    % Train L independent logistic regressions (columns)
    W = zeros(size(X_M,2), Lchips);
    for l = 1:Lchips
        y = double(Y_M(:,l));
        W(:,l) = train_logreg_lbfgs(X_M, y, 1e-3, 300);
    end

    % Predict on test set
    Pte = 1 ./ (1 + exp(-(Xte * W))); % Nte x L
    Rhat = Pte >= 0.5;
    perchip_acc(im) = mean(Rhat == Rte, 'all');

    % For each test CRP, build attack waveform and evaluate
    Nte = size(Rte,1);
    accept_cnt = 0;
    corr_accum = 0;
    for i = 1:Nte
        spread_pred = 2*double(Rhat(i,:)) - 1; % ±1
        spread_true = 2*double(Rte(i,:)) - 1;
        % correlation (normalized cosine similarity)
        c = (spread_pred * spread_true.') / (norm(spread_pred)*norm(spread_true) + eps);
        corr_accum = corr_accum + c;

        % simulate a spoofed CDMA slot with predicted spread
        tx_bits = randi([0,1], 1, cfg.bits_per_slot);
        Lloc = cfg.L;
        tx_imp = zeros(1, cfg.bits_per_slot * Lloc);
        for bb = 1:cfg.bits_per_slot
            idx = (bb-1)*Lloc + (1:Lloc);
            tx_imp(idx) = (2*tx_bits(bb)-1) * spread_pred;
        end
        y_imp = add_awgn_local(tx_imp, cfg.impostor_snr_db);
        slot_accept = decide_cdma_legit_local(y_imp, spread_true, tx_bits, Lloc, cfg.threshold, cfg.accept_ratio);
        if slot_accept, accept_cnt = accept_cnt + 1; end
    end
    attack_success(im) = accept_cnt / max(1,Nte);
    mean_corr(im) = corr_accum / max(1,Nte);
    fprintf('M=%d: per-chip=%.3f, corr=%.3f, attack=%.3f\n', M, perchip_acc(im), mean_corr(im), attack_success(im));
end

% Plot
figure('Name','Modeling: per-chip acc, corr, attack success vs M','Color','w');
yyaxis left;
plot(Mset, attack_success, '-o', 'LineWidth', 1.8, 'DisplayName','Attack success (slot)'); hold on; grid on;
ylabel('Attack success (slot)'); ylim([0 1]);
yyaxis right;
plot(Mset, perchip_acc, '--s', 'LineWidth', 1.5, 'DisplayName','Per-chip accuracy'); hold on;
plot(Mset, mean_corr, ':d', 'LineWidth', 1.4, 'DisplayName','Mean spread correlation');
ylabel('Per-chip acc / mean corr'); ylim([0 1]);
xlabel('CRP leakage M');
title(sprintf('Modeling attack metrics vs M  (L=%d, \\tau=%.2f, acc=%.2f, impostorSNR=%g dB)', cfg.L, cfg.threshold, cfg.accept_ratio, cfg.impostor_snr_db), 'Interpreter','none');
legend('Location','northwest'); legend boxoff;

% Optional CSV export
if ~isempty(cfg.save_csv)
    T = table(Mset(:), perchip_acc(:), mean_corr(:), attack_success(:), ...
        'VariableNames', {'M','PerChipAcc','MeanCorr','AttackSuccess'});
    try
        writetable(T, cfg.save_csv);
        fprintf('Saved results to %s\n', cfg.save_csv);
    catch
        warning('Failed to write CSV: %s', cfg.save_csv);
    end
end

end

% === Helpers ===
function [Cmat, Rmat] = pack_crp_logs_expanded(CRP_LOGS, nStages, L)
N = numel(CRP_LOGS); Cmat = false(N,nStages); Rmat = false(N,L);
for i = 1:N
    Ci = logical(CRP_LOGS{i}.C(:).');
    Ri = logical(CRP_LOGS{i}.R(:).');
    if numel(Ci) < nStages, Ci = [Ci, false(1, nStages-numel(Ci))]; end
    if numel(Ci) > nStages, Ci = Ci(1:nStages); end
    if numel(Ri) < L, Ri = [Ri, false(1, L-numel(Ri))]; end
    if numel(Ri) > L, Ri = Ri(1:L); end
    Cmat(i,:) = Ci; Rmat(i,:) = Ri;
end
end

function Phi = apuf_features(Cmat)
B = 1 - 2*double(Cmat);              % 0->+1, 1->-1
P = fliplr(cumprod(fliplr(B),2));    % suffix products
Phi = [ones(size(Cmat,1),1), P];
end

function w = train_logreg_lbfgs(X, y, lambda, iters)
[N,d] = size(X); y = y(:); w = zeros(d,1); alpha = 0.5;
for t = 1:iters
    z = X*w; p = 1./(1+exp(-z));
    g = X.'*(p - y)/N + lambda*w;
    w = w - alpha*g;
end
end

function [slot_accept, pass_count] = decide_cdma_legit_local(rx_vec, spread_legit, tx_bits, L, thr, accept_ratio)
Nbits = length(tx_bits); pass_count = 0;
for bb = 1:Nbits
    idx = (bb-1)*L + (1:L); seg = rx_vec(idx);
    corr_legit = sum(seg .* spread_legit);
    rx_bit = (corr_legit > 0);
    conf = min(1, abs(corr_legit)/L);
    if conf >= thr && (rx_bit == tx_bits(bb)), pass_count = pass_count + 1; end
end
need = ceil(accept_ratio * Nbits);
slot_accept = (pass_count >= need);
end

function y = add_awgn_local(x, SNRdB)
sigP = mean(abs(x).^2); noiseP = sigP / (10^(SNRdB/10));
y = x + sqrt(noiseP/2) * randn(size(x));
end

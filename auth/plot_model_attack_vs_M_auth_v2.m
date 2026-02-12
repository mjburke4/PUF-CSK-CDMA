
function plot_model_attack_vs_M_auth_v2(varargin)
% PLOT_MODEL_ATTACK_VS_M_AUTH_V2
% Improved: Modeling attack success vs CRP leakage M using run_star_tdma_auth,
% with explicit control over CRP tap ('raw' or 'expanded') and seed-expansion.
%
% Key options:
%   'crp_tap'           : 'expanded' (default) | 'raw'
%   'use_seed_expansion': false (default for modeling) | true
%   'shuffle'           : true (default) shuffle CRP order before train/test split
%   'feature'           : 'apuf' (suffix-products) | 'bits' (raw challenge bits)
%
% Usage:
%   plot_model_attack_vs_M_auth_v2('Mset',[64 128 256 512], 'L',128, 'challenge_len',64, ...
%       'train_frames', 2000, 'test_frames', 400, 'SNRdB', 8, 'impostor_snr_db', 5, ...
%       'threshold',0.80, 'accept_ratio',0.80, 'bits_per_slot',128, ...
%       'crp_tap','expanded', 'use_seed_expansion', false, 'feature','apuf', 'rng_seed', 123);
%
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
addParameter(p,'crp_tap', 'expanded');           % 'expanded' | 'raw'
addParameter(p,'use_seed_expansion', false);     % IMPORTANT: false => learnable mapping C->R
addParameter(p,'feature', 'apuf');               % 'apuf' | 'bits'
addParameter(p,'rng_seed', 123);
addParameter(p,'shuffle', true);
parse(p, varargin{:});
cfg = p.Results;

rng(cfg.rng_seed);

% ---- 1) Harvest CRPs for training/testing ----
res_tr = run_star_tdma_auth( ...
    'scheme','CDMA', ...
    'N', 1, 'K_users_per_slot', 1, ...
    'L', cfg.L, 'bits_per_slot', 1, ...
    'frames', cfg.train_frames, ...
    'SNRdB_vec', 0, ...
    'use_nonce', true, 'nonce_len', cfg.challenge_len, ...
    'replay_reuse_prob', 0.0, ...
    'auth_mode','exposedCRP', 'crp_tap', cfg.crp_tap, ...
    'use_seed_expansion', cfg.use_seed_expansion, ...
    'log_level', 0, 'rng_seed', cfg.rng_seed);

res_te = run_star_tdma_auth( ...
    'scheme','CDMA', ...
    'N', 1, 'K_users_per_slot', 1, ...
    'L', cfg.L, 'bits_per_slot', 1, ...
    'frames', cfg.test_frames, ...
    'SNRdB_vec', 0, ...
    'use_nonce', true, 'nonce_len', cfg.challenge_len, ...
    'replay_reuse_prob', 0.0, ...
    'auth_mode','exposedCRP', 'crp_tap', cfg.crp_tap, ...
    'use_seed_expansion', cfg.use_seed_expansion, ...
    'log_level', 0, 'rng_seed', cfg.rng_seed + 1);

% Pack CRPs
if strcmpi(cfg.crp_tap,'expanded')
    [Ctr, Rtr] = pack_crp_logs_expanded(res_tr.CRP_LOGS, cfg.challenge_len, cfg.L);
    [Cte, Rte] = pack_crp_logs_expanded(res_te.CRP_LOGS, cfg.challenge_len, cfg.L);
    Lchips = size(Rtr,2);
else
    [Ctr, rtr] = pack_crp_logs_raw(res_tr.CRP_LOGS, cfg.challenge_len);
    [Cte, rte] = pack_crp_logs_raw(res_te.CRP_LOGS, cfg.challenge_len);
    Lchips = 1;  % raw tap gives 1-bit labels per CRP
end

% Shuffle if requested
if cfg.shuffle
    idx_tr = randperm(size(Ctr,1));
    Ctr = Ctr(idx_tr,:);
    if strcmpi(cfg.crp_tap,'expanded'), Rtr = Rtr(idx_tr,:); else, rtr = rtr(idx_tr); end
end

% Build features
if strcmpi(cfg.feature,'apuf')
    Xtr = apuf_features(Ctr);
    Xte = apuf_features(Cte);
else
    Xtr = [ones(size(Ctr,1),1), double(Ctr)];
    Xte = [ones(size(Cte,1),1), double(Cte)];
end

% ---- 2) Train and evaluate vs M ----
Mset = cfg.Mset(:).';
acc_attack = nan(size(Mset));     % acceptance rate at BS vs M
acc_pred   = nan(size(Mset));     % mean per-chip prediction accuracy

for im = 1:numel(Mset)
    M = min(Mset(im), size(Xtr,1));
    X_M = Xtr(1:M, :);

    if strcmpi(cfg.crp_tap,'expanded')
        Y_M = Rtr(1:M, :);             % [M x L]
        % Train L independent logistic models
        W = zeros(size(X_M,2), Lchips);
        for l = 1:Lchips
            y = double(Y_M(:,l));
            W(:,l) = train_logreg_lbfgs(X_M, y, 1e-2, 250);
        end
        Pte = 1 ./ (1 + exp(-(Xte * W)));
        Rhat = Pte >= 0.5;
        acc_pred(im) = mean(Rhat(:) == Rte(:));

        % --- Attack evaluation ---
        L = size(Rte,2);
        Nte = size(Xte,1);
        accept_cnt = 0;
        for i = 1:Nte
            spread_pred = 2*double(Rhat(i,:)) - 1;
            spread_true = 2*double(Rte(i,:))  - 1;
            tx_bits = randi([0,1], 1, cfg.bits_per_slot);
            tx_imp = zeros(1, cfg.bits_per_slot * L);
            for bb = 1:cfg.bits_per_slot
                idx = (bb-1)*L + (1:L);
                tx_imp(idx) = (2*tx_bits(bb)-1) * spread_pred;
            end
            y_imp = add_awgn_local(tx_imp, cfg.impostor_snr_db);
            [slot_accept, ~] = decide_cdma_legit_local(y_imp, spread_true, tx_bits, L, cfg.threshold, cfg.accept_ratio);
            if slot_accept, accept_cnt = accept_cnt + 1; end
        end
        acc_attack(im) = accept_cnt / max(1,Nte);
    else
        % RAW tap => 1-bit labels. Train one LR to predict r from C.
        y = double(rtr(1:M));
        W = train_logreg_lbfgs(X_M, y, 1e-2, 300);
        pte = 1 ./ (1 + exp(-(Xte * W)));
        rhat = (pte >= 0.5);
        acc_pred(im) = mean(rhat(:) == rte(:));

        % For attack with RAW tap, we can only spoof a single-bit “spread”.
        % Emulate by repeating the same chip across L (weak attack but illustrative).
        L = cfg.L;
        Nte = size(Xte,1);
        accept_cnt = 0;
        for i = 1:Nte
            chip = 2*double(rhat(i)) - 1;
            spread_pred = chip * ones(1, L);
            spread_true = 2*double(rte(i)) - 1;   % this is only 1 chip; not directly comparable
            % Since true L-chip spread is unknown in RAW tap, skip strict slot test and
            % approximate acceptance by correlation sign agreement across bits.
            % (You can extend this with a mapping from raw r to code index if desired.)
            tx_bits = randi([0,1], 1, cfg.bits_per_slot);
            tx_imp = zeros(1, cfg.bits_per_slot * L);
            for bb = 1:cfg.bits_per_slot
                idx = (bb-1)*L + (1:L);
                tx_imp(idx) = (2*tx_bits(bb)-1) * spread_pred;
            end
            y_imp = add_awgn_local(tx_imp, cfg.impostor_snr_db);
            % Here we cannot compute the *true* spread; report the mean per-chip accuracy instead
            % and set attack success as NaN to avoid misleading interpretation.
            % (If you have a defined mapping from raw r -> code index/permutation, plug it here.)
            % For now:
            accept_cnt = NaN; Nte = NaN;
        end
        acc_attack(im) = NaN;
    end
end

% ---- Plot ----
figure('Name','Modeling attack vs M (AUTH v2)','Color','w'); hold on; grid on;
plot(Mset, acc_attack, '-o', 'LineWidth', 1.8, 'DisplayName', 'Attack success');
yyaxis right;
plot(Mset, acc_pred, '--s', 'LineWidth', 1.5, 'DisplayName', 'Mean per-chip accuracy');
ylabel('Per-chip accuracy'); ylim([0.5 1]);
yyaxis left; ylabel('Attack success rate'); ylim([0 1]);
xlabel('CRP leakage M');
title(sprintf('Modeling attack vs M  (SNR=%g dB, impostor SNR=%g dB, N=%d, tap=%s, seedexp=%d, \\tau=%.2f, acc=%.2f)', ...
    cfg.SNRdB, cfg.impostor_snr_db, cfg.L, cfg.crp_tap, cfg.use_seed_expansion, cfg.threshold, cfg.accept_ratio), 'Interpreter','none');
legend('Location','northwest'); legend boxoff;
end


% ===================== Helpers =====================
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

function [Cmat, rvec] = pack_crp_logs_raw(CRP_LOGS, nStages)
N = numel(CRP_LOGS); Cmat = false(N,nStages); rvec = false(N,1);
for i = 1:N
    Ci = logical(CRP_LOGS{i}.C(:).');
    if isfield(CRP_LOGS{i},'r'), ri = logical(CRP_LOGS{i}.r); else, ri = false; end
    if numel(Ci) < nStages, Ci = [Ci, false(1, nStages-numel(Ci))]; end
    if numel(Ci) > nStages, Ci = Ci(1:nStages); end
    Cmat(i,:) = Ci; rvec(i) = ri;
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

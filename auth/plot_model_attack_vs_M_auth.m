
function plot_model_attack_vs_M_auth(varargin)
% PLOT_MODEL_ATTACK_VS_M_AUTH
% Modeling attack success vs CRP leakage M, using run_star_tdma_auth in 'exposedCRP' mode.
%
% Flow:
%   1) Harvest training CRPs (challenges C and expanded spreads R) for a single target node.
%   2) Train L independent logistic models to predict each chip R(:,l) from C features.
%   3) Evaluate attack success at a chosen SNR by transmitting with predicted spread and
%      running the paper's accept/reject test at the BS (which uses the *true* spread).
%
% Usage:
%   plot_model_attack_vs_M_auth('Mset',[64 128 256 512], 'L',128, 'challenge_len',64, ...
%       'train_frames', 1200, 'test_frames', 300, 'SNRdB', 8, 'impostor_snr_db', 5, ...
%       'threshold',0.80, 'accept_ratio',0.80, 'feature','apuf', 'rng_seed', 123);
%
% Name-Value options:
%   'Mset'            : vector of training CRP budgets (default [64 128 256 512 1024])
%   'L'               : spreading length (default 128)
%   'challenge_len'   : #stages in APUF challenge (default 64)
%   'train_frames'    : total frames to harvest for training (default 2000)  (N=1 => #CRPs ≈ frames)
%   'test_frames'     : frames to harvest for testing (default 400)
%   'SNRdB'           : SNR at BS for evaluation (default 8)
%   'impostor_snr_db' : SNR of attacker link (default 5)
%   'threshold'       : decision threshold τ (default 0.80)
%   'accept_ratio'    : slot vote ratio (default 0.80)
%   'bits_per_slot'   : payload bits per slot (default 128)
%   'feature'         : 'apuf' (suffix-product features) | 'bits' (raw 0/1) (default 'apuf')
%   'rng_seed'        : RNG seed (default 123)
%
% Notes:
%   • We set N=1 (single target node) and K=1 to simplify harvesting.
%   • Requires run_star_tdma_auth.m on path and that it logs expanded spreads when
%     'auth_mode'='exposedCRP' and 'crp_tap'='expanded'.
%
% Outputs:
%   A figure showing attacker acceptance vs M, with optional training accuracy diagnostic.
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
addParameter(p,'feature', 'apuf');      % 'apuf' | 'bits'
addParameter(p,'rng_seed', 123);
parse(p, varargin{:});
cfg = p.Results;

rng(cfg.rng_seed);

% --------- 1) Harvest CRPs for training and testing (N=1, K=1) ----------
res_tr = run_star_tdma_auth( ...
    'scheme','CDMA', ...
    'N', 1, 'K_users_per_slot', 1, ...
    'L', cfg.L, 'bits_per_slot', 1, ...   % 1 bit per slot is enough for CRP logging
    'frames', cfg.train_frames, ...
    'SNRdB_vec', 0, ...
    'use_nonce', true, 'nonce_len', cfg.challenge_len, ...
    'replay_reuse_prob', 0.0, ...
    'auth_mode','exposedCRP', 'crp_tap','expanded', ...
    'log_level', 0, 'rng_seed', cfg.rng_seed);

res_te = run_star_tdma_auth( ...
    'scheme','CDMA', ...
    'N', 1, 'K_users_per_slot', 1, ...
    'L', cfg.L, 'bits_per_slot', 1, ...
    'frames', cfg.test_frames, ...
    'SNRdB_vec', 0, ...
    'use_nonce', true, 'nonce_len', cfg.challenge_len, ...
    'replay_reuse_prob', 0.0, ...
    'auth_mode','exposedCRP', 'crp_tap','expanded', ...
    'log_level', 0, 'rng_seed', cfg.rng_seed + 1);

% Pack CRPs into matrices (C: [Ncrp x nStages], R: [Ncrp x L])
[Ctr, Rtr] = pack_crp_logs_expanded(res_tr.CRP_LOGS, cfg.challenge_len, cfg.L);
[Cte, Rte] = pack_crp_logs_expanded(res_te.CRP_LOGS, cfg.challenge_len, cfg.L);

nStages = size(Ctr,2); L = size(Rtr,2);

% Features
if strcmpi(cfg.feature,'apuf')
    Xtr = apuf_features(Ctr);   % [Ntr x (nStages+1)]
    Xte = apuf_features(Cte);   % [Nte x (nStages+1)]
else
    Xtr = [ones(size(Ctr,1),1), double(Ctr)];  % bias + raw bits
    Xte = [ones(size(Cte,1),1), double(Cte)];
end

% --------- 2) Train L independent logistic models and evaluate vs M ---------
Mset = cfg.Mset(:).';
acc_attack = nan(size(Mset));    % attacker acceptance rate vs M
acc_pred   = nan(size(Mset));    % optional: mean per-chip prediction accuracy (diagnostic)

for im = 1:numel(Mset)
    M = Mset(im);
    M = min(M, size(Xtr,1));  % cap to available CRPs

    X_M = Xtr(1:M, :);
    Y_M = Rtr(1:M, :);            % 0/1 labels for each chip

    % Train per-chip models
    W = zeros(size(X_M,2), L);
    for l = 1:L
        y = double(Y_M(:,l));     % 0/1
        W(:,l) = train_logreg_lbfgs(X_M, y, 1e-2, 200);  % lambda=1e-2, iters=200
    end

    % Predict on test challenges
    Pte = 1 ./ (1 + exp(-(Xte * W)));     % probabilities
    Rhat = Pte >= 0.5;                    % predicted spread bits
    acc_pred(im) = mean(Rhat(:) == Rte(:));

    % --------- 3) Evaluate attack success at BS decision ---------
    % Build attacker waveforms with predicted spread; BS correlates with true spread.
    Nte = size(Xte,1);
    accept_cnt = 0;
    for i = 1:Nte
        spread_pred = 2*double(Rhat(i,:)) - 1;    % ±1
        spread_true = 2*double(Rte(i,:))  - 1;    % ±1

        % Random payload bits for the attack trial
        tx_bits = randi([0,1], 1, cfg.bits_per_slot);

        % Build impostor waveform
        tx_imp = zeros(1, cfg.bits_per_slot * L);
        for bb = 1:cfg.bits_per_slot
            idx = (bb-1)*L + (1:L);
            tx_imp(idx) = (2*tx_bits(bb)-1) * spread_pred;
        end
        % Add AWGN at impostor link SNR
        y_imp = add_awgn_local(tx_imp, cfg.impostor_snr_db);

        % BS decision: uses TRUE spread and known tx_bits
        [slot_accept, ~] = decide_cdma_legit_local(y_imp, spread_true, tx_bits, L, cfg.threshold, cfg.accept_ratio);
        if slot_accept, accept_cnt = accept_cnt + 1; end
    end
    acc_attack(im) = accept_cnt / max(1,Nte);
end

% --------- Plot ---------
figure('Name','Modeling attack success vs M (AUTH)','Color','w'); hold on; grid on;
plot(Mset, acc_attack, '-o', 'LineWidth', 1.8, 'DisplayName', 'Modeling attack success');
yyaxis right;
plot(Mset, acc_pred, '--s', 'LineWidth', 1.5, 'DisplayName', 'Mean per-chip prediction acc');
ylabel('Per-chip accuracy'); ylim([0.5 1]);
yyaxis left; ylabel('Attack success rate'); ylim([0 1]);
xlabel('CRP leakage M (training samples)');
title(sprintf('Modeling attack vs M (SNR=%g dB, impostor SNR=%g dB, N=%d, \\tau=%.2f, acc=%.2f)', ...
    cfg.SNRdB, cfg.impostor_snr_db, L, cfg.threshold, cfg.accept_ratio), 'Interpreter','none');
legend('Location','northwest'); legend boxoff;

end

% ===================== Helpers =====================

function [Cmat, Rmat] = pack_crp_logs_expanded(CRP_LOGS, nStages, L)
% Pack CRP logs (expanded tap) into matrices.
N = numel(CRP_LOGS);
Cmat = false(N, nStages);
Rmat = false(N, L);
for i = 1:N
    Ci = logical(CRP_LOGS{i}.C(:).');
    Ri = logical(CRP_LOGS{i}.R(:).');
    % pad/trim as needed
    if numel(Ci) < nStages, Ci = [Ci, false(1, nStages-numel(Ci))]; end
    if numel(Ci) > nStages, Ci = Ci(1:nStages); end
    if numel(Ri) < L, Ri = [Ri, false(1, L-numel(Ri))]; end
    if numel(Ri) > L, Ri = Ri(1:L); end
    Cmat(i,:) = Ci;
    Rmat(i,:) = Ri;
end
end

function Phi = apuf_features(Cmat)
% Build APUF "suffix product" features for each challenge row.
% Cmat: [N x nStages] in {0,1}
B = 1 - 2*double(Cmat);              % 0->+1, 1->-1
P = fliplr(cumprod(fliplr(B),2));    % suffix products
Phi = [ones(size(Cmat,1),1), P];     % add bias
end

function w = train_logreg_lbfgs(X, y, lambda, iters)
% Simple logistic regression trainer with L2 (uses gradient descent fallback if fminunc unavailable).
% X: [N x d], y: [N x 1] in {0,1}
[N,d] = size(X);
y = y(:);
w = zeros(d,1);
alpha = 0.5;                 % step size
for t = 1:iters
    z = X*w;
    p = 1./(1+exp(-z));
    g = X.'*(p - y)/N + lambda*w;      % gradient
    w = w - alpha*g;
end
end

function [slot_accept, pass_count] = decide_cdma_legit_local(rx_vec, spread_legit, tx_bits, L, thr, accept_ratio)
Nbits = length(tx_bits);
pass_count = 0;
for bb = 1:Nbits
    idx = (bb-1)*L + (1:L);
    seg = rx_vec(idx);
    corr_legit = sum(seg .* spread_legit);
    rx_bit = (corr_legit > 0);
    conf = min(1, abs(corr_legit)/L);
    if conf >= thr && (rx_bit == tx_bits(bb)), pass_count = pass_count + 1; end
end
need = ceil(accept_ratio * Nbits);
slot_accept = (pass_count >= need);
end

function y = add_awgn_local(x, SNRdB)
sigP = mean(abs(x).^2);
noiseP = sigP / (10^(SNRdB/10));
y = x + sqrt(noiseP/2) * randn(size(x));
end

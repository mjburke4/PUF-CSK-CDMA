function results = ablate_lr_four_cases(varargin)
% ABLATE_LR_FOUR_CASES
% Runs 4 ablations (expansion x nonce), harvests CRPs with your simulator,
% trains Logistic Regression per node, and plots a bar chart.
%
% Optional args (name/value):
%   'N'              (default 8)
%   'frames'         (default 10000)  % ≈80k total CRPs (~10k per node with N=8)
%   'bits_per_slot'  (default 1)
%   'SNRdB'          (default 0)
%   'L'              (default 128)
%   'challenge_len'  (default 64)
%   'rng_seed'       (default 7)
%   'quiet'          (default true)
%
% Returns struct with fields: table (4x), per_node_acc (cell), and plotting handles.

% ---- parse
ip = inputParser;
addParameter(ip,'N',8);
addParameter(ip,'frames',10000);
addParameter(ip,'bits_per_slot',1);
addParameter(ip,'SNRdB',0);
addParameter(ip,'L',128);
addParameter(ip,'challenge_len',64);
addParameter(ip,'rng_seed',7);
addParameter(ip,'quiet',true);
addParameter(ip,'use_synth_raw_case1', true);   % if true, case #1 uses pure APUF baseline
addParameter(ip,'synth_sigma_intra', 0.0);      % noise on synthetic delta (keep 0 for clean labels)
addParameter(ip,'synth_crps_per_node', 10000);  % CRPs per node for the synthetic baseline
parse(ip,varargin{:});
P = ip.Results;

cases = { ...
  struct('name','Raw (no nonce)','use_seed_expansion',false,'use_nonce',false), ...
  struct('name','Raw + Nonce','use_seed_expansion',false,'use_nonce',true), ...
  struct('name','Expansion only','use_seed_expansion',true,'use_nonce',false), ...
  struct('name','Expansion + Nonce','use_seed_expansion',true,'use_nonce',true) ...
};

acc_means = zeros(1, numel(cases));
acc_stds  = zeros(1, numel(cases));
per_node_acc = cell(1, numel(cases));

fprintf('== LR ablation (N=%d, frames=%d, bpslot=%d, SNR=%.1f dB) ==\n', ...
    P.N, P.frames, P.bits_per_slot, P.SNRdB);

for k = 1:numel(cases)
    if k == 1 && P.use_synth_raw_case1
        fprintf('>> Case #1 using SYNTHETIC RAW APUF (no expansion/nonce)\n');
        [Cmat, rvec, meta_s] = harvest_synth_raw_apuf( ...
            P.N, P.synth_crps_per_node, P.challenge_len, ...
            'sigma_intra', P.synth_sigma_intra, 'rng_seed', P.rng_seed);
        meta = struct('node_ids', meta_s.node_ids, 'nStages', meta_s.nStages);
    else
        fprintf('>> Case #%d using SIMULATOR (expansion=%d, nonce=%d)\n', ...
                k, C.use_seed_expansion, C.use_nonce);
        [Cmat, rvec, meta] = harvest_crps(P.N, P.frames, P.bits_per_slot, P.SNRdB, ...
            'L', P.L, ...
            'challenge_len', P.challenge_len, ...
            'use_nonce', C.use_nonce, ...
            'nonce_len', 64, ...
            'nonce_scope', 'per-user', ...
            'K_users_per_slot', 1, ...
            'replay_reuse_prob', 0.0, ...
            'rng_seed', P.rng_seed, ...
            'log_soft_metrics', false, ...
            'quiet', P.quiet, ...
            'use_seed_expansion', C.use_seed_expansion);
    end


    C = cases{k};
    if ~P.quiet
        fprintf('  -> %s  [expansion=%d, nonce=%d]\n', C.name, C.use_seed_expansion, C.use_nonce);
    end

    % Harvest CRPs for this case
    [Cmat, rvec, meta] = harvest_crps(P.N, P.frames, P.bits_per_slot, P.SNRdB, ...
        'L', P.L, ...
        'challenge_len', P.challenge_len, ...
        'use_nonce', C.use_nonce, ...
        'nonce_len', 64, ...
        'nonce_scope', 'per-user', ...
        'K_users_per_slot', 1, ...
        'replay_reuse_prob', 0.0, ...
        'rng_seed', P.rng_seed, ...
        'log_soft_metrics', false, ...
        'quiet', P.quiet);

    % Compute LR accuracy per node (best-of-four feature-map alignment)
    [acc_node, acc_mean, acc_std] = local_lr_best_of_four_phi(Cmat, rvec, meta.node_ids, meta.nStages);
    acc_means(k) = acc_mean;
    acc_stds(k)  = acc_std;
    per_node_acc{k} = acc_node;

    if ~P.quiet
        fprintf('     mean=%.3f  std=%.3f  (nodes=%d)\n', acc_mean, acc_std, numel(unique(meta.node_ids)));
    end
end

% --- Plot
figure; hold on;
hb = bar(acc_means, 'FaceAlpha', 0.85);
er = errorbar(1:numel(cases), acc_means, acc_stds, '.');
er.LineWidth = 1.5;
er.Color = [0 0 0];

xticks(1:numel(cases));
xticklabels({cases{1}.name, cases{2}.name, cases{3}.name, cases{4}.name});
xtickangle(20);
ylim([0.45 1.00]);
grid on;
ylabel('LR test accuracy');
title(sprintf('LR accuracy across nonce/expansion ablations (N=%d, frames=%d, L=%d, n=%d)', ...
    P.N, P.frames, P.L, P.challenge_len));

% --- Pack results
results = struct();
results.cases = cases;
results.acc_means = acc_means;
results.acc_stds = acc_stds;
results.per_node_acc = per_node_acc;
results.params = P;
results.figure = gcf;
results.bar = hb;
results.errorbar = er;
end

% ---------------- helpers ----------------

function [acc_node, acc_mean, acc_std] = local_lr_best_of_four_phi(C_matrix, r_vector, node_ids, nStages)
% Train LR per node; for each node choose the best of 4 APUF feature-map variants:
%  suffix/prefix  x  normal/reversed bit order
nodes = unique(node_ids);
acc_node = zeros(numel(nodes),1);

for ii = 1:numel(nodes)
    nd = nodes(ii);
    sel = (node_ids == nd);
    C = C_matrix(sel,:); y = r_vector(sel);
    if sum(sel) < 1500
        acc_node(ii) = NaN; continue;
    end

    % Build four Phi variants
    Phi_suf    = apuf_phi_suffix(C);
    Phi_pre    = apuf_phi_prefix(C);
    Phi_suf_R  = apuf_phi_suffix(fliplr(C));
    Phi_pre_R  = apuf_phi_prefix(fliplr(C));

    a = [lr_acc(Phi_suf, y), lr_acc(Phi_pre, y), lr_acc(Phi_suf_R, y), lr_acc(Phi_pre_R, y)];
    acc_node(ii) = max(a);
end

acc_mean = nanmean(acc_node);
acc_std  = nanstd(acc_node);
end

function Phi = apuf_phi_suffix(C)
% Suffix-product feature map (common in literature)
% C: NxNbits logical/0-1
[M,n] = size(C);
b = 1 - 2*double(C);
p = fliplr(cumprod(fliplr(b),2));
Phi = [ones(M,1) p];
end

function Phi = apuf_phi_prefix(C)
% Prefix-product feature map (equivalent to suffix on reversed challenge)
[M,n] = size(C);
b = 1 - 2*double(C);
p = cumprod(b,2);
Phi = [ones(M,1) p];
end

function acc = lr_acc(Phi, y)
% Logistic regression with bias disabled (Phi includes bias column)
rp = randperm(size(Phi,1));
Nt = floor(0.8*size(Phi,1));
PhiTr = double(Phi(rp(1:Nt),:)); yTr = double(y(rp(1:Nt)));
PhiTe = double(Phi(rp(Nt+1:end),:)); yTe = double(y(rp(Nt+1:end)));
mdlLR = fitclinear(PhiTr, yTr, 'Learner','logistic','FitBias',false, ...
                   'Regularization','ridge','Lambda',1e-3,'Solver','lbfgs','ClassNames',[0 1]);
yhat  = predict(mdlLR, PhiTe);
acc   = mean(yhat == yTe);
end

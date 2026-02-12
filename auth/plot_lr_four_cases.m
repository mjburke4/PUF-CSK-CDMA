function out = plot_lr_four_cases(varargin)
% PLOT_LR_FOUR_CASES
% Four LR scenarios using your sim's nodes + CRP tap:
%  1) Raw                 : crp_tap='raw',      use_seed_expansion=false, use_nonce=false
%  2) Raw + Nonce         : crp_tap='raw',      use_seed_expansion=false, use_nonce=true
%  3) Expansion only      : crp_tap='expanded', use_seed_expansion=true,  use_nonce=false
%  4) Expansion + Nonce   : crp_tap='expanded', use_seed_expansion=true,  use_nonce=true
%
% Requires: run_star_tdma_auth.m (with 'auth_mode','exposedCRP','crp_tap'),
%           harvest_crps.m (forwards 'use_seed_expansion' and 'crp_tap').

% ---- knobs (override via name/value) ----
ip = inputParser;
addParameter(ip,'N',8);
addParameter(ip,'frames',10000);       % with bits_per_slot=1 → ~N*frames CRPs
addParameter(ip,'bits_per_slot',1);
addParameter(ip,'SNRdB',0);
addParameter(ip,'L',128);
addParameter(ip,'challenge_len',64);
addParameter(ip,'rng_seed',7);
addParameter(ip,'quiet',false);
parse(ip,varargin{:});
P = ip.Results;

cases = { ...
  struct('name','Raw',                'tap','raw',      'exp',false,'nonce',false), ...
  struct('name','Raw + Nonce',        'tap','raw',      'exp',false,'nonce',true ), ...
  struct('name','Expansion',          'tap','expanded', 'exp',true, 'nonce',false), ...
  struct('name','Expansion + Nonce',  'tap','expanded', 'exp',true, 'nonce',true ) ...
};

acc_means = zeros(1,4);
acc_stds  = zeros(1,4);
per_node_acc = cell(1,4);

fprintf('\n== 4-case LR plot (N=%d, frames=%d, bps=%d, SNR=%.1f dB) ==\n', ...
  P.N, P.frames, P.bits_per_slot, P.SNRdB);

for k = 1:4
  Ck = cases{k};
  if ~P.quiet
    fprintf('>> Case #%d: %-20s [tap=%-8s exp=%d nonce=%d]\n', ...
      k, Ck.name, Ck.tap, Ck.exp, Ck.nonce);
  end

  % --- harvest CRPs for this case ---
  [Cmat, rvec, meta] = harvest_crps(P.N, P.frames, P.bits_per_slot, P.SNRdB, ...
      'L', P.L, ...
      'challenge_len', P.challenge_len, ...
      'use_nonce', Ck.nonce, ...
      'use_seed_expansion', Ck.exp, ...
      'crp_tap', Ck.tap, ...
      'nonce_len', 64, 'nonce_scope','per-user', ...
      'K_users_per_slot', 1, ...
      'replay_reuse_prob', 0.0, ...
      'rng_seed', P.rng_seed, ...
      'log_soft_metrics', false, ...
      'quiet', true);

  % --- LR per node (best-of-four APUF Φ variants) ---
  [acc_node, acc_mean, acc_std] = local_best_of_four_phi_lr(Cmat, rvec, meta.node_ids, meta.nStages);
  acc_means(k) = acc_mean; acc_stds(k) = acc_std; per_node_acc{k} = acc_node;

  if ~P.quiet
    fprintf('   mean=%.3f  std=%.3f\n', acc_mean, acc_std);
  end
end

% --- Plot ---
figure; hold on; grid on;
b = bar(acc_means, 'FaceAlpha',0.9);
er = errorbar(1:4, acc_means, acc_stds, '.k'); er.LineWidth = 1.5;
xticks(1:4);
xticklabels({cases{1}.name, cases{2}.name, cases{3}.name, cases{4}.name});
xtickangle(12);
ylim([0.45 1.00]);
ylabel('LR test accuracy');
title(sprintf('LR across nonce/expansion ablations  (L=%d, n=%d)', P.L, P.challenge_len));

out = struct('cases',{cases}, 'acc_means',acc_means, 'acc_stds',acc_stds, ...
             'per_node_acc',{per_node_acc}, 'params',P);
end

% ===== Helpers: best-of-four Φ, suffix/prefix, LR =====
function [acc_node, acc_mean, acc_std] = local_best_of_four_phi_lr(C, y, node_ids, nStages)
nodes = unique(node_ids);
acc_node = nan(numel(nodes),1);
for ii = 1:numel(nodes)
    nd = nodes(ii);
    sel = (node_ids == nd);
    if sum(sel) < 1500, acc_node(ii) = NaN; continue; end
    Cn = C(sel,:); yn = y(sel);
    Phi_s  = apuf_phi_suffix(Cn);
    Phi_p  = apuf_phi_prefix(Cn);
    Phi_sR = apuf_phi_suffix(fliplr(Cn));
    Phi_pR = apuf_phi_prefix(fliplr(Cn));
    a = [lr_acc(Phi_s,yn), lr_acc(Phi_p,yn), lr_acc(Phi_sR,yn), lr_acc(Phi_pR,yn)];
    acc_node(ii) = max(a);
end
acc_mean = nanmean(acc_node);
acc_std  = nanstd(acc_node);
end

function Phi = apuf_phi_suffix(C)
[N,n] = size(C);
b = 1 - 2*double(C);
p = fliplr(cumprod(fliplr(b),2));
Phi = [ones(N,1) p];
end

function Phi = apuf_phi_prefix(C)
[N,n] = size(C);
b = 1 - 2*double(C);
p = cumprod(b,2);
Phi = [ones(N,1) p];
end

function acc = lr_acc(Phi, y)
rp = randperm(size(Phi,1));
Nt = floor(0.8*size(Phi,1));
PhiTr = double(Phi(rp(1:Nt),:));  yTr = double(y(rp(1:Nt)));
PhiTe = double(Phi(rp(Nt+1:end),:));  yTe = double(y(rp(Nt+1:end)));
mdl = fitclinear(PhiTr, yTr, 'Learner','logistic','FitBias',false, ...
                 'Regularization','ridge','Lambda',1e-3, ...
                 'Solver','lbfgs','ClassNames',[0 1]);
yhat = predict(mdl, PhiTe);
acc  = mean(yhat == yTe);
end

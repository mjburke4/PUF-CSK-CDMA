function out = plot_lr_svm_four_cases(varargin)
% PLOT_LR_SVM_FOUR_CASES
% Four scenarios using your sim's nodes + CRP tap; plots LR and SVM side-by-side:
%  1) Raw                 : tap='raw',      use_seed_expansion=false, use_nonce=false
%  2) Raw + Nonce         : tap='raw',      use_seed_expansion=false, use_nonce=true
%  3) Expansion only      : tap='expanded', use_seed_expansion=true,  use_nonce=false
%  4) Expansion + Nonce   : tap='expanded', use_seed_expansion=true,  use_nonce=true
%
% Requires: run_star_tdma_auth.m with 'auth_mode','exposedCRP','crp_tap'
%           harvest_crps.m that forwards 'use_seed_expansion' and 'crp_tap'

% ---- knobs ----
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

meansLR = zeros(1,4); stdsLR = zeros(1,4);
meansSVM= zeros(1,4); stdsSVM= zeros(1,4);
per_node_LR  = cell(1,4);
per_node_SVM = cell(1,4);

fprintf('\n== 4-case LR/SVM plot (N=%d, frames=%d, bps=%d, SNR=%.1f dB) ==\n', ...
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

  % --- per-node: best-of-4 Φ for LR and for SVM (independently) ---
  [acc_node_LR, acc_node_SVM, mLR, sLR, mSVM, sSVM] = ...
      best_of_four_phi_lr_svm(Cmat, rvec, meta.node_ids, meta.nStages);

  meansLR(k) = mLR;  stdsLR(k) = sLR;  per_node_LR{k}  = acc_node_LR;
  meansSVM(k)= mSVM; stdsSVM(k)= sSVM; per_node_SVM{k} = acc_node_SVM;

  if ~P.quiet
    fprintf('   LR:  mean=%.3f  std=%.3f | SVM: mean=%.3f  std=%.3f\n', mLR, sLR, mSVM, sSVM);
  end
end

% --- grouped bar plot (LR vs SVM) ---
figure; hold on; grid on;
vals = [meansLR(:), meansSVM(:)];                    % 4 x 2
hb = bar(vals, 'grouped');                           % two bars per case
% error bars
ngroups = size(vals,1); nbars = size(vals,2);
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = hb(i).XEndPoints;
end
errorbar(x(1,:), meansLR,  stdsLR,  '.k', 'LineWidth',1.2);
errorbar(x(2,:), meansSVM, stdsSVM, '.k', 'LineWidth',1.2);
xticklabels({cases{1}.name, '', cases{2}.name, '', cases{3}.name, '', cases{4}.name});
xtickangle(12);
ylim([0.45 1.00]);
ylabel('Test accuracy');
legend({'LR','SVM'}, 'Location','southoutside', 'Orientation','horizontal');
title(sprintf('LR vs SVM across nonce/expansion ablations  (L=%d, n=%d)', P.L, P.challenge_len));

out = struct('cases',{cases}, ...
             'meansLR',meansLR, 'stdsLR',stdsLR, 'per_node_LR',{per_node_LR}, ...
             'meansSVM',meansSVM, 'stdsSVM',stdsSVM, 'per_node_SVM',{per_node_SVM}, ...
             'params',P);
end

% ===== Helpers =====

function [acc_node_LR, acc_node_SVM, mLR, sLR, mSVM, sSVM] = best_of_four_phi_lr_svm(C, y, node_ids, nStages)
nodes = unique(node_ids);
acc_node_LR  = nan(numel(nodes),1);
acc_node_SVM = nan(numel(nodes),1);

for ii = 1:numel(nodes)
    nd = nodes(ii);
    sel = (node_ids == nd);
    if sum(sel) < 1500, continue; end
    Cn = C(sel,:); yn = y(sel);

    Phi_s  = apuf_phi_suffix(Cn);
    Phi_p  = apuf_phi_prefix(Cn);
    Phi_sR = apuf_phi_suffix(fliplr(Cn));
    Phi_pR = apuf_phi_prefix(fliplr(Cn));

    % LR accuracies for 4 Φ variants
    aLR = [ lr_acc(Phi_s,yn),  lr_acc(Phi_p,yn), ...
            lr_acc(Phi_sR,yn), lr_acc(Phi_pR,yn) ];
    % SVM accuracies for 4 Φ variants
    aSVM= [ svm_acc(Phi_s,yn),  svm_acc(Phi_p,yn), ...
            svm_acc(Phi_sR,yn), svm_acc(Phi_pR,yn) ];

    acc_node_LR(ii)  = max(aLR);
    acc_node_SVM(ii) = max(aSVM);
end

mLR  = nanmean(acc_node_LR);   sLR  = nanstd(acc_node_LR);
mSVM = nanmean(acc_node_SVM);  sSVM = nanstd(acc_node_SVM);
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
% Logistic regression; Phi already includes bias → disable FitBias
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

function acc = svm_acc(Phi, y)
% Linear SVM with standardization; y in {0,1} → map to {-1,+1}
rp = randperm(size(Phi,1));
Nt = floor(0.8*size(Phi,1));
PhiTr = double(Phi(rp(1:Nt),:));  yTr = double(y(rp(1:Nt)));
PhiTe = double(Phi(rp(Nt+1:end),:));  yTe = double(y(rp(Nt+1:end)));

yTr_pm = 2*yTr - 1;  % {-1,+1}
mdl = fitcsvm(PhiTr, yTr_pm, 'KernelFunction','linear', 'Standardize',true);
yhat_pm = predict(mdl, PhiTe);
yhat = (yhat_pm > 0);
acc  = mean(yhat == yTe);
end

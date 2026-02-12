function out = ablate_lr_four_cases_min(varargin)
% Minimal, self-contained 4-case LR ablation:
%  (1) SYNTHETIC raw APUF (no nonce)
%  (2) SIM raw + Nonce
%  (3) SIM expansion only
%  (4) SIM expansion + Nonce

% ---------- args ----------
ip = inputParser;
addParameter(ip,'N',8);
addParameter(ip,'frames',10000);
addParameter(ip,'bits_per_slot',1);
addParameter(ip,'SNRdB',0);
addParameter(ip,'L',128);
addParameter(ip,'challenge_len',64);
addParameter(ip,'rng_seed',7);
addParameter(ip,'quiet',false);
addParameter(ip,'synth_crps_per_node',12000);
addParameter(ip,'synth_sigma_intra',0.0);
parse(ip,varargin{:});
P = ip.Results;

cases = { ...
  struct('name','Raw (SYNTHETIC)','src','synth','use_seed_exp',false,'use_nonce',false), ...
  struct('name','Raw+Nonce (SIM)','src','sim',   'use_seed_exp',false,'use_nonce',true), ...
  struct('name','Expansion (SIM)','src','sim',   'use_seed_exp',true, 'use_nonce',false), ...
  struct('name','Expand+Nonce (SIM)','src','sim','use_seed_exp',true, 'use_nonce',true) ...
};

acc_means = zeros(1,4); acc_stds = zeros(1,4);
per_node_acc = cell(1,4);

fprintf('\n== 4-case LR ablation (N=%d, frames=%d, bps=%d, SNR=%.1f dB) ==\n', ...
  P.N, P.frames, P.bits_per_slot, P.SNRdB);

cases = { ...
  struct('name','Raw (SYNTHETIC)','src','synth_raw','use_seed_exp',false,'use_nonce',false), ...
  struct('name','Raw+Nonce (SYNTHETIC)','src','synth_raw_nonce','use_seed_exp',false,'use_nonce',true), ...
  struct('name','Expansion (SIM)','src','sim','use_seed_exp',true,'use_nonce',false), ...
  struct('name','Expand+Nonce (SIM)','src','sim','use_seed_exp',true,'use_nonce',true) ...
};

for k = 1:4
  Ck = cases{k};
  fprintf('>> Case #%d: %s  [src=%s, exp=%d, nonce=%d]\n', ...
    k, Ck.name, Ck.src, Ck.use_seed_exp, Ck.use_nonce);

  switch Ck.src
    case 'synth_raw'
      % Case #1: pure raw APUF, no nonce
      [Cmat, rvec, meta] = harvest_synth_raw_apuf(P.N, P.synth_crps_per_node, P.challenge_len, ...
          'sigma_intra', P.synth_sigma_intra, 'rng_seed', P.rng_seed);
      node_ids = meta.node_ids; nStages = meta.nStages;

    case 'synth_raw_nonce'
      % Case #2: raw APUF but nonce-mixed challenge: C_eff = C_seed ⊕ PRG(U)
      [Cmat, rvec, node_ids, nStages] = harvest_synth_raw_apuf_with_nonce( ...
          P.N, P.synth_crps_per_node, P.challenge_len, ...
          'nonce_len', 64, 'sigma_intra', P.synth_sigma_intra, 'rng_seed', P.rng_seed);

    case 'sim'
      % Cases #3–4: simulator (expansion on/off already handled in your sim)
      [Cmat, rvec, metaSim] = harvest_crps(P.N, P.frames, P.bits_per_slot, P.SNRdB, ...
          'L', P.L, 'challenge_len', P.challenge_len, ...
          'use_nonce', Ck.use_nonce, 'use_seed_expansion', Ck.use_seed_exp, ...
          'nonce_len', 64, 'nonce_scope','per-user', ...
          'K_users_per_slot', 1, 'replay_reuse_prob', 0.0, ...
          'rng_seed', P.rng_seed, 'log_soft_metrics', false, 'quiet', true);
      node_ids = metaSim.node_ids; nStages = metaSim.nStages;
  end

  [acc_node, acc_mean, acc_std] = best_of_four_phi_lr(Cmat, rvec, node_ids, nStages);
  acc_means(k) = acc_mean; acc_stds(k) = acc_std; per_node_acc{k} = acc_node;
  fprintf('   mean=%.3f  std=%.3f\n', acc_mean, acc_std);
end


% ---------- Plot ----------
figure; hold on;
b = bar(acc_means, 'FaceAlpha',0.9); grid on;
er = errorbar(1:4, acc_means, acc_stds, '.k'); er.LineWidth = 1.5;
xticks(1:4); xticklabels({cases{1}.name,cases{2}.name,cases{3}.name,cases{4}.name});
xtickangle(20); ylim([0.45 1.0]); ylabel('LR test accuracy');
title(sprintf('LR across nonce/expansion ablations (N=%d, frames=%d, L=%d, n=%d)', ...
    P.N, P.frames, P.L, P.challenge_len));

out = struct('cases',{cases}, 'acc_means',acc_means, 'acc_stds',acc_stds, ...
             'per_node_acc',{per_node_acc}, 'params',P);
end

% ======== helpers (self-contained) ========

function [C_matrix, r_vector, meta] = harvest_synth_raw_apuf(N_nodes, crps_per_node, nStages, varargin)
ip = inputParser;
addParameter(ip,'sigma_intra',0.0);
addParameter(ip,'rng_seed',123);
parse(ip,varargin{:}); P = ip.Results;
rng(P.rng_seed,'twister');
N_total = N_nodes * crps_per_node;
C_matrix = false(N_total, nStages);
r_vector = false(N_total, 1);
node_ids = zeros(N_total,1);
W = randn(nStages+1, N_nodes); % per-node weights
row = 0;
for nd = 1:N_nodes
    C = randi([0 1], crps_per_node, nStages, 'uint8');
    Phi = apuf_phi_suffix(C);
    delta = Phi * W(:,nd) + P.sigma_intra*randn(crps_per_node,1);
    r = delta >= 0;
    C_matrix(row+1:row+crps_per_node,:) = logical(C);
    r_vector(row+1:row+crps_per_node)   = logical(r);
    node_ids(row+1:row+crps_per_node)   = nd;
    row = row + crps_per_node;
end
meta = struct('node_ids',node_ids,'nStages',nStages);
end

function [acc_node, acc_mean, acc_std] = best_of_four_phi_lr(C, y, node_ids, nStages)
nodes = unique(node_ids); acc_node = nan(numel(nodes),1);
for ii = 1:numel(nodes)
    nd = nodes(ii); sel = (node_ids==nd);
    if sum(sel) < 1500, continue; end
    Cn = C(sel,:); yn = y(sel);
    Phi_s = apuf_phi_suffix(Cn);
    Phi_p = apuf_phi_prefix(Cn);
    Phi_sR= apuf_phi_suffix(fliplr(Cn));
    Phi_pR= apuf_phi_prefix(fliplr(Cn));
    a = [lr_acc(Phi_s,yn), lr_acc(Phi_p,yn), lr_acc(Phi_sR,yn), lr_acc(Phi_pR,yn)];
    acc_node(ii) = max(a);
end
acc_mean = nanmean(acc_node); acc_std = nanstd(acc_node);
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
rp = randperm(size(Phi,1)); Nt = floor(0.8*size(Phi,1));
PhiTr = double(Phi(rp(1:Nt),:)); yTr = double(y(rp(1:Nt)));
PhiTe = double(Phi(rp(Nt+1:end),:)); yTe = double(y(rp(Nt+1:end)));
mdlLR = fitclinear(PhiTr, yTr, 'Learner','logistic','FitBias',false, ...
    'Regularization','ridge','Lambda',1e-3,'Solver','lbfgs','ClassNames',[0 1]);
yhat = predict(mdlLR, PhiTe);
acc = mean(yhat == yTe);
end

function [C_eff, r_vector, node_ids, nStages] = harvest_synth_raw_apuf_with_nonce( ...
        N_nodes, crps_per_node, nStages, varargin)
% Raw APUF (k=1) labels on nonce-mixed effective challenges C_eff = C_seed ⊕ PRG(U).
ip = inputParser;
addParameter(ip,'nonce_len',64);
addParameter(ip,'sigma_intra',0.0);
addParameter(ip,'rng_seed',321);
parse(ip,varargin{:}); P = ip.Results;

rng(P.rng_seed,'twister');

N_total = N_nodes * crps_per_node;
C_seed   = false(N_total, nStages);
C_eff    = false(N_total, nStages);
r_vector = false(N_total, 1);
node_ids = zeros(N_total,1);

% per-node APUF weights
W = randn(nStages+1, N_nodes);

row = 0;
for nd = 1:N_nodes
    C = randi([0 1], crps_per_node, nStages, 'uint8');  % seed challenges (enrollment)
    % per-CRP nonce → PRG mask
    U  = randi([0 1], crps_per_node, P.nonce_len, 'uint8');
    mask = prg_bits_from_nonce(U, nStages);             % [crps x nStages], logical

    Ce = xor(logical(C), mask);                         % effective challenges
    Phi = apuf_phi_suffix(Ce);
    delta = Phi * W(:,nd) + P.sigma_intra*randn(crps_per_node,1);
    r = delta >= 0;

    idx = row+1:row+crps_per_node;
    C_seed(idx,:)   = logical(C);
    C_eff(idx,:)    = logical(Ce);
    r_vector(idx)   = logical(r);
    node_ids(idx)   = nd;
    row = row + crps_per_node;
end
end

function M = prg_bits_from_nonce(Ubits, nStages)
% Simple, deterministic per-row PRG mask from a binary nonce.
% Ubits: [N x L_u] 0/1 → returns M: [N x nStages] 0/1
[N,Lu] = size(Ubits);
M = false(N, nStages);
for i=1:N
    % fold nonce bits to a 32-bit seed
    s = uint32(2166136261);
    for j=1:Lu
        s = bitxor(s, uint32(Ubits(i,j) + j));
        s = uint32(mod(uint64(s) * 16777619, 2^32));
    end
    x = s;
    for k=1:nStages
        x = bitxor(x, bitshift(x,13,'uint32'));
        x = bitxor(x, bitshift(x,-17,'uint32'));
        x = bitxor(x, bitshift(x,5,'uint32'));
        M(i,k) = bitand(x, uint32(1)) ~= 0;
        x = bitxor(x, uint32(2654435761));
    end
end
end

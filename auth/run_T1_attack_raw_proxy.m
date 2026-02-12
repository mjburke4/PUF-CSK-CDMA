function out = run_T1_attack_raw_proxy(varargin)
% RUN_T1_ATTACK_RAW_PROXY — Raw APUF baseline: learn 1-bit, forge by repeating bit across L
p = inputParser;
addParameter(p,'N',8); addParameter(p,'frames',5000);
addParameter(p,'bits_per_slot',1); addParameter(p,'SNRdB',0);
addParameter(p,'L',128); addParameter(p,'challenge_len',64);
addParameter(p,'rng_seed',17); addParameter(p,'thr',0.8); addParameter(p,'accept_ratio',0.8);
addParameter(p,'model','lr'); addParameter(p,'quiet',false);
parse(p,varargin{:}); P=p.Results; rng(P.rng_seed);

% Harvest RAW CRPs (one bit) from your sim
[Cmat, rvec, meta] = harvest_crps(P.N, P.frames, P.bits_per_slot, P.SNRdB, ...
  'L',P.L,'challenge_len',P.challenge_len, ...
  'use_nonce', false, 'use_seed_expansion', false, ...
  'crp_tap','raw','rng_seed',P.rng_seed,'quiet',true);

Phi_all = apuf_phi_suffix(Cmat);
nodes = unique(meta.node_ids); nN=numel(nodes);
per_node_acc = NaN(nN,1);

if ~P.quiet, fprintf('Training raw APUF (1-bit) models and simulating proxy forging...\n'); end
accepted=0; trials=0;

% recreate nodes for spread_legit
try, sim_nodes = make_nodes(P.N,P.challenge_len,128,256,false,0.7); catch, sim_nodes=[]; end

for ii=1:nN
  nd = nodes(ii);
  idx = find(meta.node_ids==nd); rp=randperm(numel(idx));
  Nt=floor(0.8*numel(idx)); tr=idx(rp(1:Nt)); te=idx(rp(Nt+1:end));
  PhiTr = double(Phi_all(tr,:)); yTr=double(rvec(tr));
  PhiTe = double(Phi_all(te,:)); yTe=double(rvec(te));

  switch lower(P.model)
    case 'lr'
      mdl = fitclinear(PhiTr,yTr,'Learner','logistic','FitBias',false,'Regularization','ridge','Lambda',1e-3,'Solver','lbfgs','ClassNames',[0 1]);
      yhat = predict(mdl,PhiTe);
    case 'svm'
      mdl = fitcsvm(PhiTr,2*yTr-1,'KernelFunction','linear','Standardize',true);
      yhat = predict(mdl,PhiTe) > 0;
  end
  acc = mean(yhat==yTe); per_node_acc(ii)=acc;

  % simulate forging on test set
  if ~isempty(sim_nodes), node = sim_nodes(nd); else, node = []; end
  for t=1:numel(te)
    C = Cmat(te(t),:);
    % predict on that C
    PhiC = apuf_phi_suffix(C);
    switch lower(P.model)
      case 'lr',    rhat = predict(mdl,double(PhiC));
      case 'svm',   rhat = predict(mdl,double(PhiC))>0;
    end
    spread_pred = (2*double(rhat)-1) * ones(1,P.L);           % repeat bit
    tmpl = build_cdma_template_from_spread(spread_pred, 1, P.L);
    if ~isempty(node),
     %Rleg = get_puf_bits(node,C,P.L); 
     r_true = rvec(te(t)); 
    else, Rleg = get_puf_bits([],C,P.L); end
    %spread_legit = 2*double(Rleg)-1;
    spread_legit = (2*double(r_true)-1) * ones(1, P.L);
    rx = add_awgn(tmpl,P.SNRdB);
    [ok,~]=decide_cdma_legit(rx,spread_legit,1,P.L,P.thr,P.accept_ratio);
    trials=trials+1; if ok, accepted=accepted+1; end
  end
end

out.per_node_mean_acc = per_node_acc;
out.overall_ml_acc = nanmean(per_node_acc);
out.attack_accept_rate = accepted/max(1,trials);
out.total_trials = trials;

if ~P.quiet
  fprintf('Raw baseline: mean ML acc=%.3f, attack_accept=%.3f (trials=%d)\n', out.overall_ml_acc, out.attack_accept_rate, trials);
end
end

function Phi = apuf_phi_suffix(C)
C = double(C); [N,~]=size(C); b=1-2*C; p=fliplr(cumprod(fliplr(b),2)); Phi=[ones(N,1) p];
end

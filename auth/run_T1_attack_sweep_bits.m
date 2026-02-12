function out = run_T1_attack_sweep_bits(varargin)
% Train once (expanded CRPs), then sweep bits_per_slot and plot attack acceptance.
% Usage:
%   out = run_T1_attack_sweep_bits('N',8,'frames',5000,'L',128,'SNRdB',0, ...
%                                  'challenge_len',64,'thr',0.8,'accept_ratio',0.8, ...
%                                  'bits_list',[1 8 32 64 128]);

% ---- args ----
ip = inputParser;
addParameter(ip,'N',8);
addParameter(ip,'frames',5000);
addParameter(ip,'SNRdB',0);
addParameter(ip,'L',128);
addParameter(ip,'challenge_len',64);
addParameter(ip,'rng_seed',13);
addParameter(ip,'model','lr');                 % 'lr' or 'svm'
addParameter(ip,'thr',0.8);
addParameter(ip,'accept_ratio',0.8);
addParameter(ip,'bits_list',[1 8 32 64 128]);  % sweep
addParameter(ip,'quiet',false);
parse(ip,varargin{:});
P = ip.Results; rng(P.rng_seed);

% ---- 1) Harvest expanded CRPs once (fast harvest: use bits_per_slot=1) ----
if ~P.quiet, fprintf('Harvesting expanded CRPs (train once)...\n'); end
[Cmat, Rmat, meta] = harvest_crps(P.N, P.frames, 1, P.SNRdB, ...
    'L',P.L, 'challenge_len',P.challenge_len, ...
    'use_nonce',true, 'use_seed_expansion',true, ...
    'crp_tap','expanded', 'rng_seed',P.rng_seed, 'quiet', true);
% Defensive: rebuild Rmat if harvest_crps didn’t return expanded spreads
% if isempty(Rmat) || size(Rmat,2) == 1
%     if ~isfield(meta,'node_ids')
%         error('Need meta.node_ids to rebuild expanded spreads.');
%     end
%     % Recreate node structs deterministically (must match make_nodes defaults)
%     try
%         recon_nodes = make_nodes(P.N, P.challenge_len, 128, 256, false, 0.7);
%     catch
%         error('Cannot rebuild Rmat: make_nodes not available with expected signature.');
%     end
%     S = size(Cmat,1);
%     Rmat = false(S, P.L);
%     for s = 1:S
%         nd = meta.node_ids(s);                      % 1..N
%         Rmat(s,:) = get_puf_bits(recon_nodes(nd), Cmat(s,:), P.L);
%     end
% end
node_ids = meta.node_ids(:);  nodes = unique(node_ids); nN = numel(nodes);

% Feature map (suffix-phi)
Phi_all = apuf_phi_suffix(Cmat);

% ---- 2) Train per-node, per-chip models once ----
if ~P.quiet, fprintf('Training per-chip %s models (L=%d)...\n', P.model, P.L); end
models = cell(nN, P.L); test_idx_store = cell(nN,1);
for ii = 1:nN
    nd = nodes(ii);
    idx = find(node_ids==nd);
    rp = randperm(numel(idx));
    Nt = floor(0.8*numel(idx));
    tr = idx(rp(1:Nt)); te = idx(rp(Nt+1:end));
    test_idx_store{ii} = te;

    PhiTr = double(Phi_all(tr,:));
    for j = 1:P.L
        yTr = double(Rmat(tr,j));
        if numel(unique(yTr)) < 2, models{ii,j} = []; continue; end
        switch lower(P.model)
            case 'lr'
                mdl = fitclinear(PhiTr, yTr, 'Learner','logistic','FitBias',false, ...
                                 'Regularization','ridge','Lambda',1e-3,'Solver','lbfgs','ClassNames',[0 1]);
            case 'svm'
                mdl = fitcsvm(PhiTr, yTr, 'KernelFunction','linear','Standardize',true);
        end
        models{ii,j} = mdl;
    end
end

% Prepare legit nodes (for expected spread)
try
    sim_nodes = make_nodes(P.N, P.challenge_len, 128, 256, false, 0.7);
catch
    sim_nodes = [];
end

% ---- 3) Sweep bits_per_slot ----
acc_rate = zeros(1, numel(P.bits_list));
for b = 1:numel(P.bits_list)
    B = P.bits_list(b);
    if ~P.quiet, fprintf('Evaluating forged acceptance (bits_per_slot=%d)...\n', B); end
    total_trials = 0; accepted = 0;

    for ii = 1:nN
        nd = nodes(ii); te = test_idx_store{ii};
        if isempty(te), continue; end
        if ~isempty(sim_nodes), node = sim_nodes(nd); else, node = []; end

        for t = 1:numel(te)
            s = te(t);
            C = Cmat(s,:);

            % Predict full spread for this challenge
            PhiC = apuf_phi_suffix(C);
            Rhat = zeros(1,P.L);
            for j = 1:P.L
                mdl = models{ii,j};
                if isempty(mdl), Rhat(j) = 0; continue; end
                switch lower(P.model)
                    case 'lr',  Rhat(j) = predict(mdl, double(PhiC));
                    case 'svm', Rhat(j) = double(predict(mdl, double(PhiC)) > 0);
                end
            end
            spread_pred = 2*double(Rhat) - 1;     % 1 x L

            % Build forged slot with B bits (repeat spread per bit)
            tx_bits = ones(1, B);                 % bit value doesn't matter for correlation
            tmpl_forged = build_cdma_template_from_spread(spread_pred, tx_bits, P.L);

            % Legit expected spread from true node PUF
            % if ~isempty(node), R_legit = get_puf_bits(node, C, P.L);
            % else,              R_legit = get_puf_bits([],  C, P.L);
            % end
            % spread_legit = 2*double(R_legit) - 1;

            % (NEW) use the truth from the log for this sample s
            R_true      = Rmat(s, :);             % 1 x L, 0/1
            spread_legit = 2*double(R_true) - 1;  % ±1

            % Pass through AWGN at verifier SNR
            rx = add_awgn(tmpl_forged, P.SNRdB);

            % Verifier decision with vote over B bits
            [ok, ~] = decide_cdma_legit(rx, spread_legit, tx_bits, P.L, P.thr, P.accept_ratio);
            total_trials = total_trials + 1;
            if ok, accepted = accepted + 1; end
        end
    end
    acc_rate(b) = accepted / max(1,total_trials);
end

% ---- 4) Plot + return ----
figure; plot(P.bits_list, acc_rate, '-o','LineWidth',1.8); grid on;
xlabel('bits\_per\_slot'); ylabel('Attack acceptance rate');
title(sprintf('T1 forged acceptance vs bits/slot  (L=%d, thr=%.2f, acc=%.2f)', P.L, P.thr, P.accept_ratio));

out = struct('bits_list',P.bits_list, 'accept_rate',acc_rate, 'params',P);
end

% --- helper (suffix phi) ---
function Phi = apuf_phi_suffix(C)
C = double(C);
[N,~]=size(C); b = 1-2*C; p = fliplr(cumprod(fliplr(b),2)); Phi = [ones(N,1) p];
end

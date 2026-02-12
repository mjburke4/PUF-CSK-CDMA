function out = run_T1_attack_expanded(varargin)
% RUN_T1_ATTACK_EXPANDED  Exact (PHY-level) T1 demo: harvest expanded CRPs,
% train per-chip ML (LR by default), build forged templates from predicted
% spreads, and evaluate verifier acceptance rate.
%
% Usage:
%   out = run_T1_attack_expanded('N',8,'frames',5000,'L',128,'SNRdB',0);
%
% Returns struct out with fields:
%   out.per_chip_acc      : node_count x L matrix of per-chip test accuracies
%   out.per_node_mean_acc : node_count x 1 mean accuracy across chips
%   out.overall_ml_acc    : scalar mean over all nodes and chips (test set)
%   out.attack_accept_rate: scalar (forged acceptance fraction)
%   out.total_trials      : number of forging trials
%   out.params            : struct of parameters used
%
% Notes:
%  - Requires your harvest_crps to return expanded R vectors (Rmat) when
%    called with 'use_seed_expansion',true and 'crp_tap','expanded'.
%  - Uses suffix-phi features by default. Uses per-chip logistic regression.
%  - If you want SVM instead, pass 'model','svm' (slower).
%  - For speed, reduce frames or N during experimentation.

% -------------------- Parse inputs --------------------
p = inputParser;
addParameter(p,'N',8);
addParameter(p,'frames',5000);
addParameter(p,'bits_per_slot',1);
addParameter(p,'SNRdB',0);
addParameter(p,'L',128);
addParameter(p,'challenge_len',64);
addParameter(p,'rng_seed',13);
addParameter(p,'train_frac',0.8);
addParameter(p,'model','lr'); % 'lr' or 'svm'
addParameter(p,'quiet',false);
addParameter(p,'thr',0.8);
addParameter(p,'accept_ratio',0.8);

parse(p,varargin{:});
P = p.Results;

% Pack params for return
out = struct();
out.params = P;

% basic prints
if ~P.quiet
    fprintf('RUN_T1_ATTACK_EXPANDED: N=%d frames=%d bits_per_slot=%d L=%d SNR=%.1f dB model=%s thr=%.2f acc=%.2f\n', ...
        P.N, P.frames, P.bits_per_slot, P.L, P.SNRdB, P.model, P.thr, P.accept_ratio);
end

rng(P.rng_seed);

% -------------------- 1) Harvest expanded CRPs --------------------
if ~P.quiet, fprintf('Harvesting expanded CRPs (this may take a while)...\n'); end
% Expected returns:
%   Cmat: S x nC (challenges)
%   Rmat: S x L (expanded response/spread bits for each sample)
%   meta.node_ids: S x 1 (node index per sample)
[Cmat, Rmat, meta] = harvest_crps(P.N, P.frames, P.bits_per_slot, P.SNRdB, ...
    'L', P.L, 'challenge_len', P.challenge_len, ...
    'use_nonce', true, 'use_seed_expansion', true, ...
    'crp_tap', 'expanded', 'rng_seed', P.rng_seed, 'quiet', P.quiet);

if isempty(Cmat) || isempty(Rmat)
    error('harvest_crps did not return expected Cmat/Rmat. Ensure crp_tap=''expanded'' and harvest_crps supports it.');
end

S = size(Cmat,1);
if ~isfield(meta,'node_ids')
    error('harvest_crps must return meta.node_ids (node index per sample).');
end
node_ids = meta.node_ids(:);
nodes_list = unique(node_ids);
nNodes = numel(nodes_list);
L = size(Rmat,2);
if L ~= P.L
    warning('Returned R length (%d) differs from requested L (%d). Using returned L.', L, P.L);
    P.L = L; out.params.L = L;
end

% Build feature matrix Phi (Suffix-phi). Can be changed to best-of-four outside.
Phi_all = apuf_phi_suffix(Cmat);  % S x F

% -------------------- 2) Train per-node, per-chip models --------------------
train_frac = P.train_frac;
per_chip_acc = NaN(nNodes, P.L);   % node x chip
models = cell(nNodes, P.L);
test_idx_store = cell(nNodes,1);

if ~P.quiet, fprintf('Training per-node, per-chip models (%s)...\n', P.model); end
for ni = 1:nNodes
    nd = nodes_list(ni);
    idx = find(node_ids == nd);
    Nsamp = numel(idx);
    if Nsamp < 100
        if ~P.quiet, fprintf(' node %d: only %d samples (skip)\n', nd, Nsamp); end
        continue;
    end

    rp = randperm(Nsamp);
    Nt = max(1, floor(train_frac * Nsamp));
    tr_local = idx(rp(1:Nt));
    te_local = idx(rp(Nt+1:end));
    test_idx_store{ni} = te_local;

    PhiTr = double(Phi_all(tr_local,:));
    PhiTe = double(Phi_all(te_local,:));

    % Loop chips
    for j = 1:P.L
        yTr = double(Rmat(tr_local, j));
        yTe = double(Rmat(te_local, j));

        % Skip trivial labels
        if numel(unique(yTr)) < 2
            % degenerate labels in train set: mark accuracy as NaN
            per_chip_acc(ni,j) = NaN;
            models{ni,j} = [];
            continue;
        end

        switch lower(P.model)
            case 'lr'
                mdl = fitclinear(PhiTr, yTr, 'Learner','logistic', 'FitBias', false, ...
                                 'Regularization','ridge','Lambda',1e-3, 'Solver','lbfgs','ClassNames',[0 1]);
                yhat = predict(mdl, PhiTe);
            case 'svm'
                mdl = fitcsvm(PhiTr, yTr, 'KernelFunction','linear', 'Standardize',true);
                yhat_pm = predict(mdl, PhiTe);
                yhat = double(yhat_pm > 0);
            otherwise
                error('Unknown model type: %s', P.model);
        end

        acc = mean(yhat == yTe);
        per_chip_acc(ni,j) = acc;
        models{ni,j} = mdl;
    end

    if ~P.quiet
        fprintf(' node %2d: trained on %4d samples, test %4d -- mean chip acc = %.4f\n', ...
            nd, numel(tr_local), numel(te_local), nanmean(per_chip_acc(ni,:)));
    end
end

per_node_mean_acc = nanmean(per_chip_acc,2); % mean across chips for each node
overall_ml_acc = nanmean(per_node_mean_acc);

% -------------------- 3) Simulate forging & verifier interaction --------------------
if ~P.quiet, fprintf('Simulating PHY-level forging and verifier decisions...\n'); end
total_trials = 0;
accepted = 0;

% If make_nodes exists, recreate nodes with same factory parameters so get_puf_bits works.
% We assume challenge_len and puf_bits defaults match harvest_crps usage.
try
    nodes_struct = make_nodes(P.N, P.challenge_len, 128, 256, false, 0.7);
catch
    % If make_nodes not available with same signature, set nodes_struct empty and rely on get_puf_bits implementing logic
    nodes_struct = [];
end

for ni = 1:nNodes
    nd = nodes_list(ni);
    te_idx = test_idx_store{ni};
    if isempty(te_idx), continue; end
    % If node struct exists, pick it
    if ~isempty(nodes_struct)
        node = nodes_struct(nd);
    else
        node = []; % get_puf_bits should accept node as used elsewhere
    end

    % For each test sample predict whole spread using per-chip models
    for t = 1:numel(te_idx)
        s = te_idx(t);
        C = Cmat(s,:);
        % Predict R_hat
        Rhat = zeros(1, P.L);
        for j = 1:P.L
            mdl = models{ni,j};
            if isempty(mdl)
                % fallback: predict majority class from training set (or 0)
                Rhat(j) = 0;
            else
                PhiC = apuf_phi_suffix(C);  % 1 x F
                if strcmpi(P.model,'lr')
                    Rhat(j) = predict(mdl, double(PhiC));
                else
                    ypm = predict(mdl, double(PhiC));
                    Rhat(j) = double(ypm > 0);
                end
            end
        end

        % Build forged template (exact PHY waveform)
        spread_pred = 2*double(Rhat) - 1;  % ±1
        tx_bits = 1;          % bits_per_slot was 1 in harvest; adapt if different
        tmpl_forged = build_cdma_template_from_spread(spread_pred, tx_bits, P.L); % 1 x L

        % Legitimate spread (verifier's expected)
        if ~isempty(node)
            R_legit = get_puf_bits(node, C, P.L);
        else
            % If node not available, call get_puf_bits using a placeholder node struct
            % We assume get_puf_bits can accept a node-like struct or seed. If not, user should ensure make_nodes is available.
            R_legit = get_puf_bits([], C, P.L);
        end
        spread_legit = 2*double(R_legit) - 1;

        % Transmit forged template through AWGN (verifier SNR)
        rx = add_awgn(tmpl_forged, P.SNRdB);

        % Verifier decision
        [slot_accept, ~] = decide_cdma_legit(rx, spread_legit, tx_bits, P.L, P.thr, P.accept_ratio);

        total_trials = total_trials + 1;
        if slot_accept
            accepted = accepted + 1;
        end
    end
end

attack_accept_rate = accepted / max(1, total_trials);

% -------------------- 4) Prepare outputs --------------------
out.per_chip_acc = per_chip_acc;            % nNodes x L
out.per_node_mean_acc = per_node_mean_acc;  % nNodes x 1
out.overall_ml_acc = overall_ml_acc;
out.attack_accept_rate = attack_accept_rate;
out.total_trials = total_trials;

% Print summary
if ~P.quiet
    fprintf('\n---- Summary ----\n');
    fprintf('Mean ML accuracy across nodes (mean per-node mean) = %.4f\n', out.overall_ml_acc);
    fprintf('Attack acceptance rate (forged PHY-level) = %.4f (trials=%d)\n', attack_accept_rate, total_trials);
    fprintf('Per-node mean accuracies:\n');
    for ni = 1:nNodes
        fprintf(' node %2d: mean chip acc = %.4f\n', nodes_list(ni), per_node_mean_acc(ni));
    end
    fprintf('-----------------\n');
end

end

% -------------------- Helper: apuf_phi_suffix --------------------
function Phi = apuf_phi_suffix(C)
% C: S x n (0/1)
C = double(C);
[N,~] = size(C);
b = 1 - 2*C;
p = fliplr(cumprod(fliplr(b),2));
Phi = [ones(N,1) p];
end

function demo_T1_attack_sim()
% DEMO_T1_ATTACK_SIM
% End-to-end T1 demo:
%  - harvest CRPs (attacker sniff)
%  - train per-node LR
%  - forge responses using LR and attempt impersonation (simulated)
%
% Assumes: harvest_crps exists and returns [Cmat, rvec, meta]
%          meta.node_ids gives node index per sample
%          run_star_tdma_auth helpers (build_cdma_template, add_awgn, decide_cdma_legit) are on path.

%% --------- User knobs (reduce frames if slow) ----------
N           = 8;            % number of nodes
frames      = 3000;         % reduce for speed (2000..5000 recommended)
bits_per_slot = 1;         % 1 for CRP harvesting speed
SNRdB       = 0;            % SNR for harvesting (attacker sniffs at this SNR)
L           = 128;          % spread length used in verifier decisions
challenge_len = 64;
rng_seed    = 11;
quiet       = false;
%% -------------------------------------------------------

rng(rng_seed);

fprintf('Harvesting CRPs: N=%d frames=%d (this may take a moment)...\n', N, frames);
[Cmat, rvec, meta] = harvest_crps(N, frames, bits_per_slot, SNRdB, ...
    'L', L, 'challenge_len', challenge_len, ...
    'use_nonce', true, 'use_seed_expansion', false, ...    % we harvest raw tap here
    'crp_tap', 'raw', ...
    'rng_seed', rng_seed, 'quiet', quiet);

% meta.node_ids should be provided; otherwise infer nodes round-robin
if isfield(meta,'node_ids')
    node_ids = meta.node_ids;
else
    % fallback: assume CRPs were produced in round-robin node order
    node_ids = repmat(1:N, 1, ceil(size(Cmat,1)/N));
    node_ids = node_ids(1:size(Cmat,1))';
end

% If meta contains Ceff (keeps original shape), use that; else use Cmat
if isfield(meta,'Ceff_list'), Ceff_list = meta.Ceff_list; else Ceff_list = Cmat; end

% Per-node train/test split indices
unique_nodes = unique(node_ids);
nNodes = numel(unique_nodes);

% Train a per-node LR classifier (C -> r) and evaluate accuracy
lr_models = cell(nNodes,1);
acc_node = nan(nNodes,1);
fprintf('Training per-node LR models...\n');
for ii = 1:nNodes
    nd = unique_nodes(ii);
    idx = find(node_ids == nd);
    X = double(Ceff_list(idx,:));    % challenges as rows (0/1)
    y = double(rvec(idx));           % labels 0/1
    % require reasonable data volume
    if numel(y) < 200, fprintf(' node %d: not enough samples (%d) - skipping\n', nd, numel(y)); continue; end

    % random split 80/20
    rp = randperm(numel(y));
    Nt = floor(0.8 * numel(y));
    tr = rp(1:Nt); te = rp(Nt+1:end);

    Phi = apuf_phi_suffix(X);   % choose suffix-phi as feature (you can replace with best-of-four)
    PhiTr = double(Phi(tr,:));  yTr = y(tr);
    PhiTe = double(Phi(te,:));  yTe = y(te);

    mdl = fitclinear(PhiTr, yTr, 'Learner','logistic', 'FitBias', false, ...
                     'Regularization','ridge','Lambda',1e-3,'Solver','lbfgs','ClassNames',[0 1]);
    lr_models{ii} = mdl;

    yhat = predict(mdl, PhiTe);
    acc = mean(yhat == yTe);
    acc_node(ii) = acc;
    fprintf(' node %2d: LR accuracy = %.4f  (samples=%d)\n', nd, acc, numel(y));
end

fprintf('Mean LR accuracy across nodes = %.4f (std=%.4f)\n', nanmean(acc_node), nanstd(acc_node));

%% -------- Simulate impersonation attempts using LR predictions ----------
% For each held-out test sample we will:
%  - predict r_hat from its C
%  - construct a forged CDMA template using spread_pred = (2*r_hat-1) * ones(1,L)
%  - pass the forged template through AWGN at the verifier SNR and run decide_cdma_legit
% Notes: this is a proxy. For exact PHY-level attack, harvest expanded spreads and train to predict those.

fprintf('Simulating impersonation attempts (proxy forging)...\n');
total_trials = 0;
accepted = 0;

% We'll reuse the same split as training stage: for each node, use its test indices
for ii = 1:nNodes
    nd = unique_nodes(ii);
    idx = find(node_ids == nd);
    if numel(idx) < 200, continue; end

    % recreate same split as before
    rp = randperm(numel(idx));
    Nt = floor(0.8 * numel(idx));
    te_idx_local = rp(Nt+1:end);
    te_idx = idx(te_idx_local);   % indices in global arrays

    % get node object for access to true get_puf_bits
    % NOTE: we assume make_nodes was called inside harvest_crps and meta has nodes W or same seed.
    % If you need the actual node struct, call make_nodes(...) with same rng / params.
    nodes = make_nodes(N, challenge_len, 128, 256, false, 0.7);  % must match earlier parameters for node mapping
    node = nodes(nd);

    mdl = lr_models{ii};
    if isempty(mdl), continue; end

    for t = 1:numel(te_idx)
        sidx = te_idx(t);
        C = Ceff_list(sidx,:);

        % Predict r_hat
        PhiC = apuf_phi_suffix(double(C));
        rhat = predict(mdl, double(PhiC));   % 0/1

        % Build forged template: repeated predicted-bit spread (proxy)
        spread_pred = (2*double(rhat)-1) * ones(1, L);   % ±1 vector of length L
        tx_bits = 1;  % single payload bit (we'll assume tx_bit=1)
        tmpl_forged = build_cdma_template_from_spread(spread_pred, tx_bits, L); % returns 1xL

        % Verifier's legitimate spread & decision: build expected spread from true PUF
        R_legit = get_puf_bits(node, C, L);
        spread_legit = 2*R_legit - 1;

        % Attacker transmits forged template over AWGN at verifier SNR: use SNRdB (same as harvest)
        rx = add_awgn(tmpl_forged, SNRdB);

        % Use decide_cdma_legit to simulate verifier's slot decision
        [slot_accept, pass_count] = decide_cdma_legit(rx, spread_legit, tx_bits, L, 0.8, 0.8);

        total_trials = total_trials + 1;
        if slot_accept, accepted = accepted + 1; end
    end
end

attack_accept_rate = accepted / max(1, total_trials);
fprintf('Simulated impersonation acceptance rate (proxy) = %.4f  (trials=%d)\n', attack_accept_rate, total_trials);

%% Summary
fprintf('\nSummary:\n - Mean LR node accuracy = %.4f\n - Proxy attack acceptance rate = %.4f\n', nanmean(acc_node), attack_accept_rate);
fprintf('Notes: this is a proxy simulation. For precise PHY-level impersonation, harvest expanded spreads (crp_tap=''expanded'') and train to predict spreads.\n');

end


%% ===== Helper: apuf_phi_suffix (same as used in your other scripts) =====
function Phi = apuf_phi_suffix(C)
C = double(C);
[N,~] = size(C);
b = 1 - 2*C;
p = fliplr(cumprod(fliplr(b),2));
Phi = [ones(N,1) p];
end

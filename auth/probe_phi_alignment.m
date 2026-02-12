function probe_phi_alignment(res_like)
% res_like.CRP_LOGS = meta.logs; res_like.config = meta.config;

logs = res_like.CRP_LOGS;
nStages = res_like.config.challenge_len;

% pack per-node
N = numel(logs);
C_all = false(N, nStages); y_all = false(N,1); node_all = zeros(N,1);
for i=1:N
    C_all(i,:)  = logs{i}.C(:).';
    y_all(i)    = logs{i}.r;
    node_all(i) = logs{i}.node;
end
nodes = unique(node_all);

fprintf('--- Probing feature-map alignment (nStages=%d) ---\n', nStages);
for nd = nodes.'
    idx = (node_all == nd);
    C = C_all(idx,:);  y = y_all(idx);
    if sum(idx) < 1000, fprintf('node %d: too few CRPs (%d)\n', nd, sum(idx)); continue; end

    % Build four feature variants
    Phi_suffix   = apuf_phi_suffix(C);                 % our original
    Phi_prefix   = apuf_phi_prefix(C);                 % prefix product
    Phi_suffix_R = apuf_phi_suffix(fliplr(C));         % reversed bits + suffix
    Phi_prefix_R = apuf_phi_prefix(fliplr(C));         % reversed bits + prefix

    % Evaluate LR on each (80/20 split)
    acc = @(Phi) lr_acc(Phi, y);
    a = [acc(Phi_suffix), acc(Phi_prefix), acc(Phi_suffix_R), acc(Phi_prefix_R)];
    [best,which] = max(a);
    names = ["suffix","prefix","suffix(revC)","prefix(revC)"];

    fprintf('node %2d: best LR=%.4f using %s  (all: %s)\n', nd, best, names(which), ...
        sprintf(' %.3f', a));
end
end

function Phi = apuf_phi_suffix(C)
% C: NxNbits (0/1). Phi: Nx(Nbits+1). (suffix product, common in papers)
[M,n] = size(C);
b = 1 - 2*double(C);                 % 0->+1, 1->-1
p = fliplr(cumprod(fliplr(b),2));    % p_i = prod_{j=i..n} b_j
Phi = [ones(M,1) p];
end

function Phi = apuf_phi_prefix(C)
% Prefix-product variant (equivalent to suffix on reversed bit order)
[M,n] = size(C);
b = 1 - 2*double(C);
p = cumprod(b, 2);                   % p_i = prod_{j=1..i} b_j
Phi = [ones(M,1) p];
end

function acc = lr_acc(Phi, y)
% Use logistic regression; Phi already contains bias column -> disable built-in bias
rp = randperm(size(Phi,1));
Nt = floor(0.8*size(Phi,1));
PhiTr = double(Phi(rp(1:Nt),:)); yTr = double(y(rp(1:Nt)));
PhiTe = double(Phi(rp(Nt+1:end),:)); yTe = double(y(rp(Nt+1:end)));
mdlLR = fitclinear(PhiTr, yTr, 'Learner','logistic','FitBias',false, ...
                   'Regularization','ridge','Lambda',1e-3,'Solver','lbfgs','ClassNames',[0 1]);
yhat  = predict(mdlLR, PhiTe);
acc   = mean(yhat == yTe);
end

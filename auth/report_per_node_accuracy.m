function report_per_node_accuracy(res_like)
% res_like.CRP_LOGS = meta.logs; res_like.config = meta.config;
logs = res_like.CRP_LOGS;
nStages = res_like.config.challenge_len;

% Gather arrays
N = numel(logs);
C_all = false(N, nStages);
y_all = false(N,1);
node_all = zeros(N,1);
for i=1:N
    C_all(i,:)  = logs{i}.C(:).';
    y_all(i)    = logs{i}.r;
    node_all(i) = logs{i}.node;
end

nodes = unique(node_all);
fprintf('Per-node LR/SVM accuracy (nStages=%d)\n', nStages);
for nd = nodes.'
    idx = (node_all == nd);
    C = C_all(idx,:); y = y_all(idx);
    Phi = apuf_phi(uint8(C));

    % split
    rp = randperm(sum(idx));
    Nt = floor(0.8*sum(idx));
    PhiTr = double(Phi(rp(1:Nt),:)); yTr = double(y(rp(1:Nt)));
    PhiTe = double(Phi(rp(Nt+1:end),:)); yTe = double(y(rp(Nt+1:end)));

    % LR
    mdlLR = fitclinear(PhiTr, yTr, 'Learner','logistic', ...
                       'Regularization','ridge','Lambda',1e-3, ...
                       'Solver','lbfgs','ClassNames',[0 1]);
    yhatLR = predict(mdlLR, PhiTe);
    accLR = mean(yhatLR == yTe);

    % linear SVM
    mdlSVM = fitcsvm(PhiTr, 2*yTr-1, 'KernelFunction','linear','Standardize',true);
    yhatSVM = predict(mdlSVM, PhiTe);
    accSVM = mean( (yhatSVM > 0) == (yTe > 0) );

    fprintf('  node %2d: LR=%.4f  SVM=%.4f  (CRPs=%d)\n', nd, accLR, accSVM, sum(idx));
end
end

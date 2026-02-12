function demo_train_lr_svm_from_results(res, nStages)
% Build X (Phi) and y from CRP_LOGS, then train/test LR/SVM.

logs = res.CRP_LOGS;  % cell array of structs
N = numel(logs);
Cmat = false(N, nStages);
y    = false(N, 1);

for i = 1:N
    Ci = logs{i}.C(:).';
    Cmat(i,1:nStages) = Ci(1:nStages);   % truncate if needed
    y(i) = logs{i}.r;
end

Phi = apuf_phi(uint8(Cmat));   % same Phi map used in APUF modeling

% Train/test split
idx = randperm(N);
Nt = floor(0.8*N);
PhiTr = double(Phi(idx(1:Nt),:)); yTr = double(y(idx(1:Nt)));
PhiTe = double(Phi(idx(Nt+1:end),:)); yTe = double(y(idx(Nt+1:end)));

% Logistic Regression
mdlLR = fitclinear(PhiTr, yTr, 'Learner','logistic', ...
                   'Regularization','ridge', 'Lambda',1e-3, ...
                   'Solver','lbfgs', 'ClassNames',[0 1]);
yhatLR = predict(mdlLR, PhiTe);
accLR  = mean(yhatLR == yTe);

% Linear SVM
mdlSVM = fitcsvm(PhiTr, 2*yTr-1, 'KernelFunction','linear', 'Standardize',true);
yhatSVM = predict(mdlSVM, PhiTe);
accSVM  = mean( (yhatSVM > 0) == (yTe > 0) );

fprintf('LR acc=%.4f, SVM acc=%.4f (N=%d, nStages=%d)\n', accLR, accSVM, N, nStages);
end

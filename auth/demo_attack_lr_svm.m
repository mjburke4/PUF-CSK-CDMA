function demo_attack_lr_svm()
% Generate CRPs from an APUF/XOR-APUF and train LR/SVM to model it.

nStages = 64; kXOR = 1; Ntrain = 20000; Ntest = 5000; sigmaIntra = 0.00;
rng(2,'twister');

% Enrollment weights
W = randn(nStages+1, kXOR);

% Training set
Ctr  = randi([0 1], Ntrain, nStages, 'uint8');
PhiT = apuf_phi(Ctr);
yT   = false(Ntrain,1);
for i=1:Ntrain
    yT(i) = xor_apuf_single(Ctr(i,:), W, kXOR, sigmaIntra);
end

% Test set
Cte  = randi([0 1], Ntest, nStages, 'uint8');
PhiE = apuf_phi(Cte);
yE   = false(Ntest,1);
for i=1:Ntest
    yE(i) = xor_apuf_single(Cte(i,:), W, kXOR, sigmaIntra);
end

% --- Logistic Regression ---
fprintf('\nLogistic Regression (fitclinear)...\n');
mdlLR = fitclinear(double(PhiT), double(yT), ...
    'Learner','logistic', 'Regularization','ridge', ...
    'Lambda',1e-3, 'Solver','lbfgs', 'ClassNames',[0 1]);
yHatLR = predict(mdlLR, double(PhiE));
accLR = mean(yHatLR == double(yE));
fprintf('  Test accuracy: %.4f\n', accLR);

% --- Linear SVM ---
fprintf('\nLinear SVM (fitcsvm)...\n');
yTpm = double(yT)*2 - 1; % map to {-1,+1}
mdlSVM = fitcsvm(double(PhiT), yTpm, ...
    'KernelFunction','linear', 'BoxConstraint',1, 'Standardize',true);
yHatSVM = predict(mdlSVM, double(PhiE));
accSVM = mean( (yHatSVM > 0) == (yTpm(1:numel(yHatSVM))>0) );
fprintf('  Test accuracy: %.4f\n', accSVM);
end

function bit = xor_apuf_single(C, W, k, sigma)
Phi = apuf_phi(C);
b = false(1,k);
for i=1:k
    delta = Phi * W(:,i) + sigma*randn(1,1);
    b(i) = (delta>=0);
end
bit = b(1);
for j=2:k, bit = xor(bit,b(j)); end
end

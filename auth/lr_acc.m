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
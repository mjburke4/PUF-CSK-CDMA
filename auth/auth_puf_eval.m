function [rBits, deltaMat] = auth_puf_eval(C, W, kXOR, sigma)
% Evaluate XOR-APUF on a single challenge C (1 x n), with noise sigma.
% Returns response bits and raw per-APUF deltas.
if isrow(C); C = C(:)'; end
Phi = apuf_phi(C);       % 1 x (n+1)
k   = size(W,2);
deltaMat = zeros(1,k);
bMat = false(1,k);
for i=1:k
    delta = Phi * W(:,i) + sigma * randn(1,1);
    deltaMat(i) = delta;
    bMat(i) = (delta >= 0);
end
% XOR the k APUF outputs
rBits = xor_reduce(bMat);
end

function y = xor_reduce(B)
% B: 1 x k logical
y = B(1);
for j=2:numel(B)
    y = xor(y,B(j));
end
end

function Phi = apuf_phi_prefix(C)
% Prefix-product feature map (equivalent to suffix on reversed challenge)
[M,n] = size(C);
b = 1 - 2*double(C);
p = cumprod(b,2);
Phi = [ones(M,1) p];
end
function Phi = apuf_phi_suffix(C)
% APUF suffix-product feature map (standard in the literature).
% C: [N x n] logical/0-1, Phi: [N x (n+1)] with leading 1 (bias).
[N,n] = size(C);
b = 1 - 2*double(C);                 % 0->+1, 1->-1
p = fliplr(cumprod(fliplr(b),2));    % p_i = Π_{j=i..n} b_j
Phi = [ones(N,1), p];
end

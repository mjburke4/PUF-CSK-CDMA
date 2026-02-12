function Phi = apuf_phi(C)
% Feature vector for additive-delay APUF model.
% C: M x n bits (0/1). Phi: M x (n+1). Phi(:,1)=1.
[M,n] = size(C);
b = 1 - 2*double(C);                 % map 0->+1, 1->-1
p = fliplr(cumprod(fliplr(b),2));    % p_i = prod_{j=i}^n b_j
Phi = [ones(M,1), p];
end

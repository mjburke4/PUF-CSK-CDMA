function C_set = derive_challenges(C_seed, K)
% derive_challenges - Expand a single seed challenge into K per-bit challenges (row-wise)
% Inputs:
%   C_seed : 1 x Lc binary seed challenge (row or column)
%   K      : number of one-bit challenges to derive
% Output:
%   C_set  : K x Lc matrix, each row a derived challenge
%
% Deterministic and symmetric: both ends derive identical challenge sets from the same seed.

    C_seed = double(C_seed(:)');     % ensure row vector
    Lc = numel(C_seed);
    % Simple reproducible scalar seed from bits (portable across MATLAB versions)
    prng_seed = sum((1:Lc) .* C_seed);
    rng(prng_seed, 'twister');
    C_set = randi([0,1], K, Lc);
end
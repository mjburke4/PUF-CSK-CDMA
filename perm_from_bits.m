function perm = perm_from_bits(R, M)
% perm_from_bits - Create a permutation of length M from PUF bits
% Inputs:
%   R : vector of PUF bits (0/1), at least ceil(log2(M!)) bits ideally
%   M : size of the permutation
%
% Output:
%   perm : permutation indices of 1:M

    % Ensure R is binary row vector
    R = R(:)';  
    if isempty(R)
        error('R must be non-empty');
    end

    % Convert chunks of bits into pseudo-random ranks
    needed_bits = ceil(log2(M * M)); % a few more bits to spread values
    if length(R) < needed_bits
        R = repmat(R, 1, ceil(needed_bits/length(R))); 
    end

    % Hash the bitstring into a deterministic numeric stream
    seed = sum(R .* (2.^(0:length(R)-1)));
    rng(seed, 'twister');

    % Generate a permutation of 1:M
    perm = randperm(M);

    % Reset RNG (so as not to affect other MATLAB randomness)
    rng('shuffle');
end

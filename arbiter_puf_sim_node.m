function R = arbiter_puf_sim_node(nd, C_vec, Nbits)
    % Require a scalar struct with field A
    if ~isstruct(nd) || ~isfield(nd,'A') || ~isscalar(nd)
        error('arbiter_puf_sim_node: first arg must be a SCALAR struct with field A (got size %s).', mat2str(size(nd)));
    end

    A = nd.A;                                   % challenge_len x puf_bits
    [challenge_len, puf_bits] = size(A);

    % Normalize challenge to +/-1 row
    if ischar(C_vec) || isstring(C_vec)
        Cb = double(char(C_vec)) - 48;          % '0'/'1' -> 0/1
    else
        Cb = double(C_vec);
    end
    Cb = Cb(:).';
    if numel(Cb) ~= challenge_len
        error('arbiter_puf_sim_node: challenge length %d != %d.', numel(Cb), challenge_len);
    end
    Cpm = Cb; Cpm(Cpm==0) = -1;

    % Project once to get puf_bits and take sign
    s = (A.' * Cpm.').';                        % 1 x puf_bits
    Rfull = (s > 0);                            % logical

    % Return exactly Nbits
    if Nbits <= puf_bits
        R = double(Rfull(1:Nbits));
    else
        reps = ceil(Nbits / puf_bits);
        R = repmat(double(Rfull), 1, reps);
        R = R(1:Nbits);
    end
end


function nodes = make_nodes(N, challenge_len, puf_bits, crps_per_node, use_node_similarity, similarity_alpha)
% make_nodes - Create N nodes with unique Arbiter PUF matrices and CRP tables
% Optional similarity model:
%   If use_node_similarity=true, each node's A is composed as:
%   A_i = alpha * G + sqrt(1 - alpha^2) * E_i
%   where G is a global shared component and E_i is node-unique noise.
%
% Args:
%   N, challenge_len, puf_bits, crps_per_node: as before
%   use_node_similarity (bool): default false
%   similarity_alpha (0..1): default 0.0 (no similarity)

    if nargin < 6, similarity_alpha = 0.0; end
    if nargin < 5, use_node_similarity = false; end
    if nargin < 4, crps_per_node = 128; end
    if nargin < 3, puf_bits = 128; end
    if nargin < 2, challenge_len = 64; end

    nodes = struct('id', {}, 'A', {}, 'CRP_C', {}, 'CRP_R', {});

    rng(12345); % deterministic build for reproducibility

    if use_node_similarity
        rng(54321);
        G = randn(challenge_len, puf_bits); % shared component
    end

    for i = 1:N
        node.id = i;
        if use_node_similarity
            rng(100000 + i);
            Ei = randn(challenge_len, puf_bits);
            a = max(0,min(1,similarity_alpha));
            node.A = a * G + sqrt(max(0,1 - a^2)) * Ei;
            % After node.A is formed (either random or with similarity)
            node.A = node.A ./ vecnorm(node.A);   % column-wise L2 normalization
        else
            node.A  = randn(challenge_len, puf_bits);  % secret matrix for this node
            % After node.A is formed (either random or with similarity)
            node.A = node.A ./ vecnorm(node.A);   % column-wise L2 normalization
        end

        % Factory enrollment: pre-generate CRPs
        node.CRP_C = randi([0,1], crps_per_node, challenge_len);
        node.CRP_R = zeros(crps_per_node, puf_bits);
        for r = 1:crps_per_node
            C = node.CRP_C(r, :);
            node.CRP_R(r, :) = arbiter_puf_sim_node(node, C, puf_bits);
        end

        nodes(i) = node;
    end
end
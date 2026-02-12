function [C_matrix, r_vector, meta] = harvest_synth_raw_apuf(N_nodes, crps_per_node, nStages, varargin)
% HARVEST_SYNTH_RAW_APUF
% Generate a *pure* APUF (k=1) dataset per node using the classical
% additive-delay model: delta = Phi(C) * w + noise; r = 1[delta>=0].
%
% Inputs:
%   N_nodes        : number of independent devices
%   crps_per_node  : CRPs per device (e.g., 10000)
%   nStages        : challenge length (e.g., 64)
% Optional (name/value):
%   'sigma_intra'  : std dev of additive Gaussian noise on delta (default 0)
%   'rng_seed'     : RNG seed (default 123)
%
% Outputs:
%   C_matrix : [N_total x nStages] logical challenges
%   r_vector : [N_total x 1]       logical responses
%   meta     : struct with node_ids, nStages, sigma_intra, etc.

ip = inputParser;
addParameter(ip,'sigma_intra',0.0);
addParameter(ip,'rng_seed',123);
parse(ip,varargin{:});
P = ip.Results;

rng(P.rng_seed,'twister');

N_total = N_nodes * crps_per_node;
C_matrix = false(N_total, nStages);
r_vector = false(N_total, 1);
node_ids = zeros(N_total,1);

% Draw one independent weight vector per node (n+1)
W = randn(nStages+1, N_nodes);  % columns = nodes

row = 0;
for nd = 1:N_nodes
    % Challenges for this node
    C = randi([0 1], crps_per_node, nStages, 'uint8');
    Phi = apuf_phi_suffix(C);                 % [crps_per_node x (n+1)]
    delta = Phi * W(:,nd) + P.sigma_intra*randn(crps_per_node,1);
    r = delta >= 0;

    C_matrix(row+1:row+crps_per_node, :) = logical(C);
    r_vector(row+1:row+crps_per_node)    = logical(r);
    node_ids(row+1:row+crps_per_node)    = nd;

    row = row + crps_per_node;
end

meta = struct();
meta.N_nodes       = N_nodes;
meta.crps_per_node = crps_per_node;
meta.nStages       = nStages;
meta.sigma_intra   = P.sigma_intra;
meta.rng_seed      = P.rng_seed;
meta.node_ids      = node_ids;
end

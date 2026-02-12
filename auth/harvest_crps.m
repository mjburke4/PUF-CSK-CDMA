function [C_matrix, r_vector, meta] = harvest_crps(N, frames, bits_per_slot, SNRdB, varargin)
% HARVEST_CRPS  One-liner CRP collector for LR/SVM modeling.
% Returns:
%   if crp_tap='expanded':
%       C_matrix : [Ncrp x nStages]
%       R_matrix : [Ncrp x L]          % <-- 2nd output is full spread
%       meta     : struct
%   if crp_tap='raw':
%       C_matrix : [Ncrp x nStages]
%       r_vector : [Ncrp x 1]
%       meta     : struct
%
% Example:
%   [C,r,meta] = harvest_crps(8, 2000, 1, 0);
%   demo_train_lr_svm_from_results(struct('CRP_LOGS', meta.logs, 'config', meta.config), meta.nStages);

    % ---------- Defaults / overrides ----------
    ip = inputParser;
    addParameter(ip,'L',128);
    addParameter(ip,'challenge_len',64);
    addParameter(ip,'use_nonce',true);
    addParameter(ip,'nonce_len',64);
    addParameter(ip,'nonce_scope','per-user');    % 'per-user' | 'broadcast'
    addParameter(ip,'rng_seed', 9001);
    addParameter(ip,'log_soft_metrics', false);   % set true to also log correlator mags
    addParameter(ip,'K_users_per_slot', 1);
    addParameter(ip,'replay_reuse_prob', 0.0);    % keep 0 for clean freshness
    addParameter(ip,'quiet', true);               % suppress sim printing
    addParameter(ip,'save_mat', '');              % e.g., 'crps_run01.mat'
    addParameter(ip,'save_csv', '');              % e.g., 'crps_run01.csv'  (C||r)
    addParameter(ip,'use_seed_expansion', true);  % pass-through to run_star_tdma_auth
    addParameter(ip,'crp_tap','expanded');   % 'raw' | 'expanded' (default)

    
    parse(ip, varargin{:});
    P = ip.Results;

    % ---------- Run the sim in CRP logging mode ----------
    res = run_star_tdma_auth( ...
        'scheme', 'CDMA', ...
        'N', N, ...
        'L', P.L, ...
        'bits_per_slot', bits_per_slot, ...
        'frames', frames, ...
        'SNRdB_vec', SNRdB, ...
        'use_nonce', P.use_nonce, ...
        'nonce_len', P.nonce_len, ...
        'nonce_scope', P.nonce_scope, ...
        'replay_reuse_prob', P.replay_reuse_prob, ...
        'K_users_per_slot', P.K_users_per_slot, ...
        'auth_mode', 'exposedCRP', ...     % <-- enable CRP logging taps
        'log_soft_metrics', P.log_soft_metrics, ...
        'log_level', double(~P.quiet), ...
        'rng_seed', P.rng_seed, ...
        'crp_tap', P.crp_tap);
    
    L = P.L;  % expected spread length
    logs = res.CRP_LOGS;                 % cell array of structs {node, C, r, snr_dB}
    nStages = res.config.challenge_len;  % from your sim
    Ncrp = numel(logs);

    % ---------- Pack into matrices ----------

logs    = res.CRP_LOGS;
nStages = res.config.challenge_len;
Ncrp    = numel(logs);
L       = P.L;  % expected spread length

% Pre-allocate common fields
C_matrix = false(Ncrp, nStages);
snr_vec  = zeros(Ncrp,1);
node_vec = zeros(Ncrp,1);

% Mode-specific allocation
isExpanded = strcmpi(P.crp_tap,'expanded');
if isExpanded
    R_matrix = false(Ncrp, L);    % full spread per sample
else
    r_vector = false(Ncrp, 1);    % single bit per sample
end

for i = 1:Ncrp
    % Challenge
    Ci = logs{i}.C(:).';
    if numel(Ci) ~= nStages
        Ci = Ci(1:min(end,nStages));
        if numel(Ci) < nStages
            Ci = [Ci, false(1, nStages-numel(Ci))];
        end
    end
    C_matrix(i,:) = logical(Ci);

    % Common meta
    snr_vec(i)  = logs{i}.snr_dB;
    node_vec(i) = logs{i}.node;

    % Payload
    if isExpanded
        if ~isfield(logs{i}, 'R')
            error('CRP_LOGS(%d) missing .R for expanded tap. Ensure run_star_tdma_auth logs .R.', i);
        end
        Ri = logs{i}.R(:).';
        if numel(Ri) ~= L
            Ri = Ri(1:min(end,L));
            if numel(Ri) < L
                Ri = [Ri, false(1, L-numel(Ri))];
            end
        end
        R_matrix(i,:) = logical(Ri);
    else
        r_vector(i) = logical(logs{i}.r);
    end
end

% --- Satisfy function outputs ---
if isExpanded
    % Caller binds 2nd output to Rmat; give them the full matrix via r_vector
    r_vector = R_matrix;           % <-- important: assigns output #2
end

% ---- Metadata ----
meta = struct();
meta.Ncrp          = Ncrp;
meta.nStages       = nStages;
meta.N_nodes       = N;
meta.bits_per_slot = bits_per_slot;
meta.frames        = frames;
meta.SNRdB         = SNRdB;
meta.config        = res.config;
meta.logs          = logs;
meta.snr_dB        = snr_vec;
meta.node_ids      = node_vec;

% ---- Optional saves ----
if ~isempty(P.save_mat)
    if isExpanded
        save(P.save_mat, 'C_matrix', 'R_matrix', 'meta', '-v7');
    else
        save(P.save_mat, 'C_matrix', 'r_vector', 'meta', '-v7');
    end
end
if ~isempty(P.save_csv)
    if isExpanded
        tbl = [double(C_matrix), double(R_matrix)];
    else
        tbl = [double(C_matrix), double(r_vector)];
    end
    writematrix(tbl, P.save_csv);
end

if ~P.quiet
    fprintf('[harvest_crps] Collected %d CRPs (nStages=%d, L=%d) at SNR=%.1f dB. N=%d, frames=%d, bpslot=%d, tap=%s\n', ...
        Ncrp, nStages, L, SNRdB, N, frames, bits_per_slot, P.crp_tap);
end

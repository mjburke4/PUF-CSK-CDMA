
function thr_star = calibrate_threshold_fixed_far_runstar(varargin)
% CALIBRATE_THRESHOLD_FIXED_FAR_RUNSTAR  Pick threshold to hit a target impostor FAR at a given SNR.
%
% Name-Value:
%   'K'                (2)
%   'L'                (128)
#   'codebook'         ('walsh'|'puf-random'|'walsh-orth')
%   'SNRdB'            (5)
%   'FARtarget'        (1e-3)
%   'thr_grid'         (0.60:0.02:0.98)
%   'frames'           (8)
%   'bits_per_slot'    (128)
%   'use_sic'          (false)
%   'use_equal_power'  (true)
%   'rng_seed'         (1337)
%
% We enable compute_impostor_far and search over 'thr_grid' to find the smallest thr
% whose impostor FAR <= FARtarget (greedy). Returns last grid value if none meet target.

p = inputParser;
addParameter(p,'K',2);
addParameter(p,'L',128);
addParameter(p,'codebook','walsh-orth');
addParameter(p,'SNRdB',5);
addParameter(p,'FARtarget',1e-3);
addParameter(p,'thr_grid',0.60:0.02:0.98);
addParameter(p,'frames',8);
addParameter(p,'bits_per_slot',128);
addParameter(p,'use_sic',false);
addParameter(p,'use_equal_power',true);
addParameter(p,'rng_seed',1337);
parse(p,varargin{:});
cfg = p.Results;

thr_star = cfg.thr_grid(end);
for thr = cfg.thr_grid
    res = run_star_tdma( ...
        'scheme','CDMA', ...
        'N', max(8, 2*cfg.K), ...
        'K_users_per_slot', cfg.K, ...
        'cdma_codebook', cfg.codebook, ...
        'use_slot_dither', true, ...
        'use_equal_power', cfg.use_equal_power, ...
        'use_sic', cfg.use_sic, ...
        'L', cfg.L, ...
        'bits_per_slot', cfg.bits_per_slot, ...
        'frames', cfg.frames, ...
        'SNRdB_vec', cfg.SNRdB, ...
        'threshold', thr, ...
        'accept_ratio', 0.8, ...
        'use_nonce', true, ...
        'replay_reuse_prob', 0.0, ...
        'compute_replay_far', false, ...
        'compute_impostor_far', true, ...
        'log_level', 0, ...
        'rng_seed', cfg.rng_seed);

    FAR_imp = res.FAR_impostor(1);
    if FAR_imp <= cfg.FARtarget
        thr_star = thr;
        break;
    end
end
end

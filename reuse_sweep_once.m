function reuse_sweep_once(scheme, reuse_probs, fixed_SNR)
% reuse_sweep_once - Sweep replay_reuse_prob at a fixed SNR and plot replay metrics.
% Usage:
%   reuse_sweep_once('CSK');                 % default reuse_probs & SNR=0 dB
%   reuse_sweep_once('CDMA', 0:0.2:1, -5);   % custom sweep and SNR

    if nargin < 1 || isempty(scheme),      scheme = 'CSK'; end
    if nargin < 2 || isempty(reuse_probs), reuse_probs = [0 0.2 0.4 0.6 0.8 1.0]; end
    if nargin < 3 || isempty(fixed_SNR),   fixed_SNR = 0; end

    % Fixed config (tweak as needed)
    N               = 8;
    M               = 16;     % used for CSK
    L               = 128;    % spreading length
    bits_per_slot   = 128;
    frames          = 50;
    SNRdB_vec       = fixed_SNR;   % single SNR point
    threshold       = 0.8;
    accept_ratio    = 0.85;
    compute_impostor_far = true;
    compute_replay_far   = true;
    use_seed_expansion   = true;
    use_node_variation   = false;
    use_node_similarity  = false;
    crps_per_node        = 500;

    ReplayFAR_cond = zeros(size(reuse_probs));
    ReplayRisk_over= zeros(size(reuse_probs));
    ImpFAR         = zeros(size(reuse_probs));

    for i = 1:numel(reuse_probs)
        rp = reuse_probs(i);
        fprintf('[Sweep] %s: reuse_prob=%.2f at SNR=%g dB\n', scheme, rp, fixed_SNR);

        res = run_star_tdma( ...
            'scheme', scheme, ...
            'N', N, ...
            'M', M, ...
            'L', L, ...
            'bits_per_slot', bits_per_slot, ...
            'frames', frames, ...
            'SNRdB_vec', SNRdB_vec, ...
            'threshold', threshold, ...
            'accept_ratio', accept_ratio, ...
            'compute_impostor_far', compute_impostor_far, ...
            'compute_replay_far',   compute_replay_far, ...
            'replay_reuse_prob',    rp, ...
            'use_seed_expansion',   use_seed_expansion, ...
            'use_node_variation',   use_node_variation, ...
            'use_node_similarity',  use_node_similarity, ...
            'crps_per_node',        crps_per_node ...
        );

        % Single SNR → scalars at index 1
        ReplayFAR_cond(i) = res.FAR_replay(1);   % conditional replay FAR
        ReplayRisk_over(i)= res.ReplayOverall(1);% overall replay compromise rate
        ImpFAR(i)         = res.FAR_impostor(1); % baseline impostor FAR
    end

    figure;
    plot(reuse_probs, ReplayFAR_cond, '-o','LineWidth',1.6); hold on; grid on;
    plot(reuse_probs, ReplayRisk_over, '-s','LineWidth',1.6);
    plot(reuse_probs, ImpFAR,         '--^','LineWidth',1.2);
    xlabel('Replay reuse probability'); ylabel('Rate (0..1)');
    legend('Replay FAR (conditional)','Overall Replay Risk','True Impostor FAR','Location','best');
    title(sprintf('Attack Metrics vs Replay Reuse Probability @ SNR=%g dB (%s)', fixed_SNR, upper(scheme)));
end

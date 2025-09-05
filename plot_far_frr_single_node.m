function plot_far_frr_single_node(scheme, snr_vec, accept_ratios, opts)
% plot_far_frr_single_node  Quick FRR/FAR vs SNR for ONE node, 3 accept ratios.
% Usage:
%   plot_far_frr_single_node('CDMA', -20:5:5, [0.6 0.75 0.9]);
%   plot_far_frr_single_node('CSK',  -20:5:5, [0.6 0.75 0.9], struct('L',128,'M',16,'frames',300));
%
% Inputs:
%   scheme         : 'CDMA' or 'CSK'
%   snr_vec        : vector of SNR(dB) to sweep (e.g., -20:5:5)
%   accept_ratios  : 1x3 vector of accept_ratio settings (e.g., [0.6 0.75 0.9])
%   opts (optional): struct of overrides:
%       .L (default 128)         chips per bit (CDMA) / spreading length (CSK)
%       .M (default 16)          CSK symbol count (power of 2; only used for CSK)
%       .frames (default 300)    Monte Carlo frames per SNR
%       .bits_per_slot (128)     payload bits per slot
%       .threshold (0.7)         per-bit/symbol confidence threshold
%       .use_sic (false)         ignored for single-user; kept for consistency
%
% Requires: run_star_tdma.m on MATLAB path.

    if nargin < 1 || isempty(scheme),       scheme = 'CDMA'; end
    if nargin < 2 || isempty(snr_vec),      snr_vec = -20:5:5; end
    if nargin < 3 || isempty(accept_ratios),accept_ratios = [0.6 0.75 0.9]; end
    if nargin < 4, opts = struct(); end

    % Defaults
    if ~isfield(opts,'L'),             opts.L = 128; end
    if ~isfield(opts,'M'),             opts.M = 16;  end
    if ~isfield(opts,'frames'),        opts.frames = 300; end
    if ~isfield(opts,'bits_per_slot'), opts.bits_per_slot = 128; end
    if ~isfield(opts,'threshold'),     opts.threshold = 0.7; end
    if ~isfield(opts,'use_sic'),       opts.use_sic = false; end

    % Common base args for a single node, single user per slot
    baseArgs = { ...
        'scheme', scheme, ...
        'N', 1, ...
        'L', opts.L, ...
        'M', opts.M, ...
        'bits_per_slot', opts.bits_per_slot, ...
        'frames', opts.frames, ...
        'SNRdB_vec', snr_vec, ...
        'threshold', opts.threshold, ...
        'use_node_variation', false, ...
        'K_users_per_slot', 1, ...
        'use_equal_power', true, ...
        'use_sic', opts.use_sic, ...
        'compute_impostor_far', false, ... % we only want legit FAR/FRR here
        'compute_replay_far',   false, ...
        'replay_reuse_prob',    0 ...
    };

    % Run three times with different accept ratios
    R = cell(1, numel(accept_ratios));
    for i = 1:numel(accept_ratios)
        ar = accept_ratios(i);
        fprintf('Running %s: accept\\_ratio=%.2f ...\n', upper(scheme), ar);
        R{i} = run_star_tdma(baseArgs{:}, 'accept_ratio', ar);
    end

    snr = R{1}.config.SNRdB_vec;

    % Plot FAR
    figure; hold on; grid on;
    colors = lines(numel(accept_ratios));
    for i = 1:numel(accept_ratios)
        plot(snr, R{i}.FAR, '-o', 'LineWidth', 1.6, 'Color', colors(i,:), ...
            'DisplayName', sprintf('accept\\_ratio=%.2f', accept_ratios(i)));
    end
    xlabel('SNR (dB)'); ylabel('FAR (0..1)');
    title(sprintf('%s: FAR vs SNR (N=1, L=%d, bits=%d, thr=%.2f)', ...
        upper(scheme), opts.L, opts.bits_per_slot, opts.threshold));
    legend('Location','best');

    % Plot FRR
    figure; hold on; grid on;
    for i = 1:numel(accept_ratios)
        plot(snr, R{i}.FRR, '-s', 'LineWidth', 1.6, 'Color', colors(i,:), ...
            'DisplayName', sprintf('accept\\_ratio=%.2f', accept_ratios(i)));
    end
    xlabel('SNR (dB)'); ylabel('FRR (0..1)');
    title(sprintf('%s: FRR vs SNR (N=1, L=%d, bits=%d, thr=%.2f)', ...
        upper(scheme), opts.L, opts.bits_per_slot, opts.threshold));
    legend('Location','best');
end

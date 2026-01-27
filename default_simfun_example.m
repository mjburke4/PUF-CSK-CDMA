
function [T_legit, T_imp] = default_simfun_example(opts)
% DEFAULT_SIMFUN_EXAMPLE  Toy simulator for T statistics (for wiring only).
% Replace with your real simulator that matches your codebase.
%
% Inputs (opts):
%   .SNRdB, .K, .N, .codebook ('walsh'|'random'), .binding ('index'|'code'), .trials
%
% Output:
%   T_legit : [opts.trials x 1] correlation statistics for legitimate transmissions
%   T_imp   : [opts.trials x 1] correlation statistics for impostor attempts
%
% This simple model assumes T ~ Normal(mu, sigma^2) where:
%   mu_legit increases with SNR and N, and decreases mildly with K;
%   mu_imp is near 0 and slightly increases for random codebooks and larger K;
%   sigma decreases with SNR and increases with K.
%
% Use this only to verify plotting/plumbing; hook to your PHY for real results.

SNRlin = 10.^(opts.SNRdB/10);
K = opts.K; N = opts.N;
isRandom = strcmpi(opts.codebook,'random');

% Mean and std models (toy)
mu_legit = 0.6 * tanh(0.35*log10(1+SNRlin)) * (1 + 0.15*log10(N/64)) * (1 - 0.08*(K-1));
mu_legit = max(mu_legit, 0); % clamp
mu_imp   = 0.02 * (1 + 0.15*isRandom) * (1 + 0.12*(K-1));

sigma_legit = 0.25 / sqrt(1+SNRlin) * (1 + 0.10*(K-1));
sigma_imp   = 0.35 / sqrt(1+SNRlin) * (1 + 0.15*(K-1));

% Binding tweaks (toy): PUF-Code gives slight separation gains at same SNR
if strcmpi(opts.binding,'code')
    mu_legit = mu_legit * 1.05;
    mu_imp   = mu_imp * 0.98;
end

T_legit = mu_legit + sigma_legit * randn(opts.trials,1);
T_imp   = mu_imp   + sigma_imp   * randn(opts.trials,1);
end

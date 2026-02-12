function [authCfg, authState, channelCfg] = init_auth_cdma()
% Configuration for CDMA+PUF authentication.

% --- Authentication configuration ---
authCfg.mode         = 'controlledPUF';  % 'controlledPUF' or 'exposedCRP'
authCfg.nStages      = 128;              % APUF stages
authCfg.kXOR         = 1;                % XOR-APUF degree (1 = plain APUF)
authCfg.sigmaIntra   = 0.00;             % PUF intra-chip noise (delay std dev)
authCfg.spreadLen    = 127;              % CDMA spreading length (chips)
authCfg.EsN0dB       = 10;               % per-chip Es/N0 for channel
authCfg.repeats      = 1;                % repeats of same handshake (0/1 recommended)
authCfg.tau          = 0.4;              % accept threshold on normalized correlation (rho > tau)

% --- Enrollment / device state (PUF weights) ---
rng(1,'twister'); % fixed seed for reproducibility
authState.W = randn(authCfg.nStages+1, authCfg.kXOR); % (n+1) x k XOR-APUF weights

% --- Channel configuration (extend to Rayleigh/multipath later) ---
channelCfg.type    = "AWGN";
channelCfg.EsN0dB  = authCfg.EsN0dB;
end

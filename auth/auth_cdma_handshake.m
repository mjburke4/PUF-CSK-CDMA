function [ok, log] = auth_cdma_handshake(cfg, state, chDefault, nodeID, varargin)
% Optional name/value pairs:
%  'challenge' : 1 x nStages (uint8/logical bits)
%  'nonce'     : uint32
%  'channel'   : struct (overrides chDefault for this call)

p = inputParser;
addParameter(p,'challenge',[]);
addParameter(p,'nonce',[]);
addParameter(p,'channel',[]);
parse(p,varargin{:});
C_in   = p.Results.challenge;
nonce_in = p.Results.nonce;
ch     = p.Results.channel;
if isempty(ch), ch = chDefault; end

n  = cfg.nStages;   L = cfg.spreadLen;
ok = false;
log = struct; log.nodeID = nodeID; log.obs = []; log.CRPs = [];

% 1) Challenge + nonce: use provided if present, else generate
if isempty(C_in)
    C = randi([0 1], 1, n, 'uint8');
else
    C = uint8(C_in(:))';  % normalize shape/type
end
if isempty(nonce_in)
    nonce = uint32(randi([0 2^32-1],1));
else
    nonce = uint32(nonce_in);
end
log.C = C; log.nonce = nonce;

% 2) Device-side PUF eval (noisy as configured)
[rbits_dev, ~] = auth_puf_eval(C, state.W, cfg.kXOR, cfg.sigmaIntra);

% Baseline exposed-CRP logging (ONLY if cfg.mode says so)
if strcmpi(cfg.mode,'exposedCRP')
    log.CRPs.C = C;
    log.CRPs.r = rbits_dev;
end

% 3) PRF-like mapping to spreading code (device)
code_dev = auth_prf_to_spread(rbits_dev, nonce, L);

% 4) Transmit over your channel model (prefer your fields if present)
% We support either EsN0dB or sigma2 in your channel struct:
if isfield(ch,'sigma2')
    sigma = sqrt(ch.sigma2);
    y = ( +1 * double(code_dev(:)) ) + sigma * randn(L,1);
elseif isfield(ch,'EsN0dB')
    y = auth_cdma_tx(+1, code_dev, ch.EsN0dB);
else
    % Fallback: use default EsN0dB in cfg/chDefault
    y = auth_cdma_tx(+1, code_dev, chDefault.EsN0dB);
end

% 5) Hub-side reproduction (noiseless PUF model)
[rbits_hub, ~] = auth_puf_eval(C, state.W, cfg.kXOR, 0.0);
code_hub = auth_prf_to_spread(rbits_hub, nonce, L);

% 6) Correlate and decide
[corrVal, symbolHat] = auth_cdma_rx(y, code_hub);
rho = corrVal / L;
log.obs.corr = corrVal; log.obs.rho = rho; log.obs.symHat = symbolHat;

ok = (rho > cfg.tau);
log.accept = ok;
end


function y = auth_cdma_tx(symbol, code, EsN0dB)
% BPSK over DS-SS: y = (symbol * code) + AWGN
L = numel(code);
chips = symbol * double(code(:));  % ±1
sigma = sqrt(0.5 * 10.^(-EsN0dB/10)); % per-chip noise std (unit chip energy)
y = chips + sigma * randn(L,1);
end

function [corrVal, symbolHat] = auth_cdma_rx(y, code)
% Correlate with local code and make a hard decision on the symbol
corrVal   = sum(y(:) .* double(code(:)));
symbolHat = sign(corrVal); if symbolHat==0, symbolHat=+1; end
end

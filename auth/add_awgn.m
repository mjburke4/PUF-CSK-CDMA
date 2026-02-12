function y = add_awgn(x, SNRdB)
    sigP = mean(abs(x).^2);
    noiseP = sigP / (10^(SNRdB/10));
    y = x + sqrt(noiseP/2) * randn(size(x));
end
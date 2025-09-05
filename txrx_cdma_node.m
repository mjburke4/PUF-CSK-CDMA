function [rx_bits, FAR, FRR, success, crps_used, obs] = txrx_cdma_node(node, C_seed, tx_bits, L, SNRdB, threshold, use_seed_expansion, accept_ratio)
    if nargin < 6, threshold = 0.8; end
    if nargin < 7, use_seed_expansion = true; end
    if nargin < 8, accept_ratio = 0.8; end
    N = length(tx_bits);
    if use_seed_expansion
        R = get_puf_bits(node, C_seed, L);  crps_used = L;
    else
        R = arbiter_puf_sim_node(node, C_seed, L); crps_used = 1;
    end
    spread = 2*R - 1;
    tx_symbols = zeros(1, N*L);
    for i = 1:N
        bit = 2*tx_bits(i) - 1;
        tx_symbols((i-1)*L+1:i*L) = bit * spread;
    end
    sigP = mean(abs(tx_symbols).^2);
    noiseP = sigP / (10^(SNRdB/10));
    noise = sqrt(noiseP/2) * randn(1, length(tx_symbols));
    rx = tx_symbols + noise;
    rx_bits = zeros(1, N);
    pass_count = 0;
    for i = 1:N
        idx = (i-1)*L + (1:L);
        seg = rx(idx);
        corr_pos = sum(seg .* spread);
        rx_bits(i) = (corr_pos > 0);
        conf = min(1, abs(corr_pos) / L);
        if conf >= threshold && (rx_bits(i) == tx_bits(i))
            pass_count = pass_count + 1;
        end
    end
    need = ceil(accept_ratio * N);
    slot_accept = (pass_count >= need);
    success = all(rx_bits == tx_bits);
    FAR = 0; FRR = 0;
    if slot_accept && ~success, FAR = 1; end
    if ~slot_accept &&  success, FRR = 1; end
    obs = struct();
    obs.rx = rx; obs.L = L; obs.Nbits = N; obs.tx_bits = tx_bits; obs.accept_ratio = accept_ratio;
end

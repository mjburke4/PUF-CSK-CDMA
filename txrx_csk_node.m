function [rx_bits, FAR, FRR, success, crps_used, obs] = txrx_csk_node(node, C_seed, tx_bits, M, L, SNRdB, threshold, use_seed_expansion, accept_ratio)
    if nargin < 7, threshold = 0.8; end
    if nargin < 8, use_seed_expansion = true; end
    if nargin < 9, accept_ratio = 0.8; end
    if mod(log2(M),1) ~= 0, error('M must be power of 2'); end
    if mod(log2(L),1) ~= 0, error('L must be power of 2'); end
    if M > L, error('M must be <= L'); end
    bps = log2(M);
    if mod(length(tx_bits), bps) ~= 0, error('Length(tx_bits) must be multiple of log2(M)'); end
    Ns = length(tx_bits)/bps;
    H = hadamard(L); codebook = H(1:M, :) ./ sqrt(L);
    if use_seed_expansion
        R = get_puf_bits(node, C_seed, M); crps_used = M;
    else
        R = arbiter_puf_sim_node(node, C_seed, M); crps_used = 1;
    end
    perm_idx = perm_from_bits(R, M);
    P = codebook(perm_idx, :);
    tx_bits_mat = reshape(tx_bits, bps, [])'; sym_idx = bi2de(tx_bits_mat, 'left-msb') + 1;
    tx_syms = P(sym_idx, :);
    tx_signal = reshape(tx_syms.', 1, []);
    sigP = mean(abs(tx_signal).^2);
    noiseP = sigP / (10^(SNRdB/10));
    noise = sqrt(noiseP/2) * randn(size(tx_signal));
    rx_signal = tx_signal + noise;
    obs.rx = rx_signal;   % 1 x (Ns*L)
    rx_bits = zeros(1, Ns*bps);
    pass_count = 0;
    for s = 1:Ns
        idx = (s-1)*L + (1:L);
        seg = rx_signal(idx);
        corrs = P * seg.';
        [~, k_hat] = max(corrs);
        correct_k = sym_idx(s);
        corr_correct = corrs(correct_k);
        maxc = max(corrs);
        conf = corr_correct / max(maxc, eps);
        if conf >= threshold && (k_hat == correct_k)
            pass_count = pass_count + 1;
        end
        bits_hat = de2bi(k_hat-1, bps, 'left-msb');
        rx_bits(1,(s-1)*bps+1:s*bps) = bits_hat;
    end
    need = ceil(accept_ratio * Ns);
    slot_accept = (pass_count >= need);
    success = all(rx_bits == tx_bits);
    FAR = 0; FRR = 0;
    if slot_accept && ~success, FAR = 1; end
    if ~slot_accept &&  success, FRR = 1; end
    obs = struct();
    %obs.rx_signal = rx_signal; obs.L=L; obs.M=M; obs.Ns=Ns; obs.bps=bps;
    obs.rx = rx_signal; obs.L=L; obs.M=M; obs.Ns=Ns; obs.bps=bps;
    obs.tx_bits=tx_bits;
    obs.sym_idx_true = sym_idx; obs.accept_ratio = accept_ratio;
end

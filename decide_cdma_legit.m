function [slot_accept, pass_count] = decide_cdma_legit(rx_vec, spread_legit, tx_bits, L, thr, accept_ratio)
% Decide accept/reject for a CDMA slot using matched-filter confidence + vote.
% rx_vec:   1 x (L*Nb) received chips for this user (residual or raw slot)
% spread_legit: 1 x L (+1/-1) legit spreading sequence for this user
% tx_bits:  1 x Nb ground-truth bits for this user
% L:        chips per bit
% thr:      confidence threshold in [0,1]
% accept_ratio: fraction of bits that must pass

    Nb = numel(tx_bits);
    pass_count = 0;
    for bb = 1:Nb
        idx = (bb-1)*L + (1:L);
        seg = rx_vec(idx);
        corr = sum(seg .* spread_legit);
        rx_hat = (corr > 0);
        conf = min(1, abs(corr)/L);
        if conf >= thr && (rx_hat == tx_bits(bb))
            pass_count = pass_count + 1;
        end
    end
    need = ceil(accept_ratio * Nb);
    slot_accept = (pass_count >= need);
end

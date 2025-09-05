% decide_csk_legit.m
function [slot_accept, pass_count] = decide_csk_legit(rx_vec, P_legit, sym_idx_true, L, thr, accept_ratio)
% rx_vec: 1 x (Ns*L) chips
% P_legit: M x L (permuted codebook for this node+challenge)
% sym_idx_true: 1 x Ns (1..M)
% L: chips/symbol
% thr: confidence threshold
% accept_ratio: fraction of symbols that must pass
    Ns = numel(sym_idx_true);
    pass_count = 0;
    for s = 1:Ns
        idx = (s-1)*L + (1:L);
        seg = rx_vec(idx);
        corrs = P_legit * seg.';
        [maxc, khat] = max(corrs);
        corr_true = corrs(sym_idx_true(s));
        conf = corr_true / max(maxc, eps);
        if (khat == sym_idx_true(s)) && (conf >= thr)
            pass_count = pass_count + 1;
        end
    end
    need = ceil(accept_ratio * Ns);
    slot_accept = (pass_count >= need);
end

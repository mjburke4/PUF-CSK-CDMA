function tmpl = build_cdma_template_from_spread(spread, tx_bits, L)
    N = numel(tx_bits);
    tmpl = zeros(1, N*L);
    for i = 1:N
        bit = 2*tx_bits(i) - 1;
        tmpl((i-1)*L + (1:L)) = bit * spread;
    end
end
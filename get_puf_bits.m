function R = get_puf_bits(node, C_vec, Nbits)
% R = get_puf_bits(node, C_vec, Nbits)
% Thin wrapper so the rest of the code can keep calling one place.
    R = arbiter_puf_sim_node(node, C_vec, Nbits);
end

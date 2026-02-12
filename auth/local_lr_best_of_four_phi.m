function [acc_node, acc_mean, acc_std] = local_lr_best_of_four_phi(C_matrix, r_vector, node_ids, nStages)
% Train LR per node; for each node choose the best of 4 APUF feature-map variants:
%  suffix/prefix  x  normal/reversed bit order
nodes = unique(node_ids);
acc_node = zeros(numel(nodes),1);

for ii = 1:numel(nodes)
    nd = nodes(ii);
    sel = (node_ids == nd);
    C = C_matrix(sel,:); y = r_vector(sel);
    if sum(sel) < 1500
        acc_node(ii) = NaN; continue;
    end

    % Build four Phi variants
    Phi_suf    = apuf_phi_suffix(C);
    Phi_pre    = apuf_phi_prefix(C);
    Phi_suf_R  = apuf_phi_suffix(fliplr(C));
    Phi_pre_R  = apuf_phi_prefix(fliplr(C));

    a = [lr_acc(Phi_suf, y), lr_acc(Phi_pre, y), lr_acc(Phi_suf_R, y), lr_acc(Phi_pre_R, y)];
    acc_node(ii) = max(a);
end

acc_mean = nanmean(acc_node);
acc_std  = nanstd(acc_node);
end
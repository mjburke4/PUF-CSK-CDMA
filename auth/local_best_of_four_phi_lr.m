function [acc_node, acc_mean, acc_std] = local_best_of_four_phi_lr(C, y, node_ids, nStages)
% Train LR per node; choose best of 4 APUF feature variants:
%   suffix/prefix  x  normal/reversed bit order
nodes = unique(node_ids);
acc_node = nan(numel(nodes),1);
for ii = 1:numel(nodes)
    nd = nodes(ii);
    sel = (node_ids == nd);
    if sum(sel) < 1500, acc_node(ii) = NaN; continue; end
    Cn = C(sel,:); yn = y(sel);
    Phi_s  = apuf_phi_suffix(Cn);
    Phi_p  = apuf_phi_prefix(Cn);
    Phi_sR = apuf_phi_suffix(fliplr(Cn));
    Phi_pR = apuf_phi_prefix(fliplr(Cn));
    a = [lr_acc(Phi_s,yn), lr_acc(Phi_p,yn), lr_acc(Phi_sR,yn), lr_acc(Phi_pR,yn)];
    acc_node(ii) = max(a);
end
acc_mean = nanmean(acc_node);
acc_std  = nanstd(acc_node);
end
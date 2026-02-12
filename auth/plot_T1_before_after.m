function out = plot_T1_before_after(varargin)
% Compare Raw (proxy) vs Expanded+Nonce (exact) — LR by default
p = inputParser;
addParameter(p,'N',8); addParameter(p,'frames',5000);
addParameter(p,'bits_per_slot',1); addParameter(p,'SNRdB',0);
addParameter(p,'L',128); addParameter(p,'challenge_len',64);
addParameter(p,'rng_seed',13); addParameter(p,'thr',0.8); addParameter(p,'accept_ratio',0.8);
addParameter(p,'model','lr'); addParameter(p,'quiet',false);
parse(p,varargin{:}); P=p.Results;

% BEFORE: Raw baseline (proxy)
raw = run_T1_attack_raw_proxy('N',P.N,'frames',P.frames,'bits_per_slot',P.bits_per_slot, ...
        'SNRdB',P.SNRdB,'L',P.L,'challenge_len',P.challenge_len, ...
        'rng_seed',P.rng_seed,'thr',P.thr,'accept_ratio',P.accept_ratio, ...
        'model',P.model,'quiet',P.quiet);

% AFTER: Expanded+Nonce exact PHY
aft = run_T1_attack_expanded('N',P.N,'frames',P.frames,'bits_per_slot',P.bits_per_slot, ...
        'SNRdB',P.SNRdB,'L',P.L,'challenge_len',P.challenge_len, ...
        'rng_seed',P.rng_seed,'model',P.model,'thr',P.thr,'accept_ratio',P.accept_ratio, ...
        'quiet',P.quiet);

% Plot
figure; yyaxis left; grid on; hold on;
bar(1, raw.overall_ml_acc, 'FaceAlpha',0.9); 
bar(2, aft.overall_ml_acc, 'FaceAlpha',0.9);
ylabel('Mean ML accuracy');
xticks([1 2]); xticklabels({'Raw (baseline)','Expanded + Nonce'});
yyaxis right;
plot([1 2], [raw.attack_accept_rate aft.attack_accept_rate], '-o','LineWidth',1.8);
ylabel('Attack acceptance rate');
title(sprintf('T1 Modeling: Before vs After  (L=%d, n=%d, thr=%.2f, acc=%.2f)', P.L, P.challenge_len, P.thr, P.accept_ratio));
legend({'ML accuracy','Attack accept rate'},'Location','southoutside');
out.raw = raw; out.after = aft; out.params=P;
end

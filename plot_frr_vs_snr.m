
function plot_frr_vs_snr(varargin)
% PLOT_FRR_VS_SNR  Helper to generate "Plot 1 — FRR vs SNR" at fixed FAR target.
% 
% Usage:
%   plot_frr_vs_snr();  % uses defaults and the default_simfun_example()
%
%   % or supply options as name-value pairs:
%   plot_frr_vs_snr('SNRdB', -5:2:20, 'Kset', [1 2 4], 'N', 128, ...
%       'Codebooks', {'walsh','random'}, 'Bindings', {'index','code'}, ...
%       'FARtarget', 1e-3, 'Trials', 2e4, 'Layout','grid', ...
%       'SimFun', @default_simfun_example);
%
% Inputs (name-value):
%   SNRdB      : vector of SNR points (dB). Default: -5:2:20
%   Kset       : vector of concurrent users. Default: [1 2 4]
%   N          : chip length. Default: 128
%   Codebooks  : {'walsh','random'} subset. Default: both
%   Bindings   : {'index','code'} subset. Default: both
%   FARtarget  : target FAR (e.g., 1e-3). Default: 1e-3
%   Trials     : number of Monte Carlo trials per point. Default: 2e4
%   Layout     : 'grid' (2x2 subplots for (binding x codebook)) or 'combined'. Default: 'grid'
%   SimFun     : handle to your simulator:
%                [T_legit, T_imp] = SimFun(opts) where opts has fields:
%                .SNRdB, .K, .N, .codebook, .binding, .trials
%
% Notes:
%  - This function DOES NOT compute PHY directly; it relies on SimFun to
%    return arrays of correlation statistics T for legitimate and impostor cases.
%  - Use helper_find_threshold.m to map impostor scores -> tau for a fixed FAR target.
%
% Output:
%  - A figure with FRR vs SNR curves for K ∈ Kset, per (binding, codebook).
%
% Author: ChatGPT helper

p = inputParser;
addParameter(p, 'SNRdB', -5:2:20);
addParameter(p, 'Kset', [1 2 4]);
addParameter(p, 'N', 128);
addParameter(p, 'Codebooks', {'walsh','random'});
addParameter(p, 'Bindings', {'index','code'});
addParameter(p, 'FARtarget', 1e-3);
addParameter(p, 'Trials', 2e4);
addParameter(p, 'Layout', 'grid'); % or 'combined'
addParameter(p, 'SimFun', @default_simfun_example);
parse(p, varargin{:});
opt = p.Results;

% build all panels
nb = numel(opt.Bindings);
nc = numel(opt.Codebooks);

if strcmpi(opt.Layout,'grid')
    figure('Name','FRR vs SNR at fixed FAR','Color','w');
end

colors = lines(numel(opt.Kset));

panel = 0;
for ib = 1:nb
    binding = opt.Bindings{ib};
    for ic = 1:nc
        codebook = opt.Codebooks{ic};
        panel = panel + 1;
        if strcmpi(opt.Layout,'grid')
            subplot(nb, nc, panel);
        else
            if panel==1, figure('Name','FRR vs SNR at fixed FAR','Color','w'); end
            hold on;
        end
        hold on;
        for ik = 1:numel(opt.Kset)
            K = opt.Kset(ik);
            FRR = nan(size(opt.SNRdB));
            for is = 1:numel(opt.SNRdB)
                SNRdB = opt.SNRdB(is);
                simopts.SNRdB = SNRdB;
                simopts.K = K;
                simopts.N = opt.N;
                simopts.codebook = codebook;
                simopts.binding = binding;
                simopts.trials = opt.Trials;

                [T_legit, T_imp] = opt.SimFun(simopts);

                % Compute tau for fixed FAR target, then FRR
                tau = helper_find_threshold(T_imp, opt.FARtarget);
                FRR(is) = mean(T_legit < tau);
            end
            plot(opt.SNRdB, FRR, '-o', 'DisplayName', sprintf('K=%d', K), 'Color', colors(ik,:));
        end

        grid on;
        xlabel('SNR (dB)');
        ylabel(sprintf('FRR at FAR=%.1e', opt.FARtarget));
        ttl = sprintf('Binding: %s | Codebook: %s | N=%d', binding_label(binding), codebook_label(codebook), opt.N);
        title(ttl, 'Interpreter','none');
        legend('Location','southwest'); legend boxoff;
        ylim([0 1]);
    end
end

if ~strcmpi(opt.Layout,'grid')
    xlabel('SNR (dB)'); ylabel(sprintf('FRR at FAR=%.1e', opt.FARtarget));
    title(sprintf('FRR vs SNR | N=%d', opt.N));
    legend('Location','southwest'); legend boxoff;
    grid on; ylim([0 1]);
end

end

function tau = helper_find_threshold(T_imp, FARtarget)
% Choose tau such that fraction of impostor scores >= tau equals FARtarget.
% We sort impostor scores ascending and pick the (1 - FARtarget) quantile.
if isempty(T_imp)
    error('T_imp is empty from simulator.');
end
T_imp = sort(T_imp(:),'ascend');
idx = max(1, min(numel(T_imp), ceil((1 - FARtarget) * numel(T_imp))));
tau = T_imp(idx);
end

function s = binding_label(binding)
switch lower(binding)
    case 'index', s = 'PUF-Index';
    case 'code',  s = 'PUF-Code';
    otherwise,    s = binding;
end
end

function s = codebook_label(codebook)
switch lower(codebook)
    case 'walsh',  s = 'Walsh-orth';
    case 'random', s = 'Random';
    otherwise,     s = codebook;
end
end

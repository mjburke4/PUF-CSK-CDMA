function out = plot_cspuf_lorawan_tradeoff(cfg)
% PLOT_CSPUF_LORAWAN_TRADEOFF
% Create plots of throughput improvement vs challenge reuse K.
% Model includes: Downlink challenge time (amortized by K and/or broadcast),
% Uplink response time, Uplink payload time, and per-slot guard/switch time.
%
% OUTPUT:
%   out.results            -> table of (mode, N, K, goodput_bps, improvement_pct)
%   out.fig_per_device     -> figure handle (Per-device reuse)
%   out.fig_broadcast      -> figure handle (Broadcast once per frame)
%   out.base_goodput_bps   -> baseline goodput (per-device, K=1)
%   out.cfg                -> resolved config (with defaults)
%
% DEFAULTS (chosen to match the plots shared in our discussion):
%   Rd = 1760 bps (LoRaWAN DR3 approx), Ru = 1760 bps
%   Lc = 256 bits (challenge), Lr = 128 bits (response), Lp = 512 bits (payload)
%   Tg = 5e-3 s guard, no extra headers/preambles/FEC overhead
%   K_values = [1 2 5 10 20 50]; N_values = [4 8 16]
%
% EXAMPLE:
%   out = plot_cspuf_lorawan_tradeoff();               % use defaults
%   cfg = struct('Rd',1760,'Ru',1760,'K_values',[1 2 5 10 20 50],'N_values',[8],'Lp',512);
%   out = plot_cspuf_lorawan_tradeoff(cfg);            % custom sweep
%
% NOTES:
%   "Per-device reuse" = a device reuses the same challenge for K of its slots.
%   "Broadcast once per frame" = one challenge covers all N devices per frame;
%   with reuse factor K, you broadcast once every K frames.

if nargin < 1, cfg = struct(); end

% ---- Resolve defaults ----------------------------------------------------
cfg = set_default(cfg, 'Rd',            1760.0);     % bps (LoRaWAN DR3 ~1.76 kbps)
cfg = set_default(cfg, 'Ru',            1760.0);     % bps
cfg = set_default(cfg, 'Lc',            256);        % bits (challenge)
cfg = set_default(cfg, 'Lr',            128);        % bits (response)
cfg = set_default(cfg, 'Lp',            512);        % bits (payload per slot)
cfg = set_default(cfg, 'Tg',            5e-3);       % seconds (guard/switch)
cfg = set_default(cfg, 'hdr_dl_bits',   0);          % extra DL header bits (if any)
cfg = set_default(cfg, 'hdr_ul_bits',   0);          % extra UL header bits (if any)
cfg = set_default(cfg, 'Tpreamble_dl',  0.0);        % DL preamble time (s)
cfg = set_default(cfg, 'Tpreamble_ul',  0.0);        % UL preamble time (s)
cfg = set_default(cfg, 'fec_dl_rate',   1.0);        % code rate (e.g., 0.8 adds 25% time)
cfg = set_default(cfg, 'fec_ul_rate',   1.0);
cfg = set_default(cfg, 'K_values',      [1 2 5 10 20 50]);
cfg = set_default(cfg, 'N_values',      [4 8 16]);
cfg = set_default(cfg, 'title_prefix',  'LoRaWAN DR3 (~1.76 kbps)');
cfg = set_default(cfg, 'verbose',       true);
cfg = set_default(cfg, 'save_figs',     false);
cfg = set_default(cfg, 'save_csv',      false);
cfg = set_default(cfg, 'save_dir',      pwd);

K_values = cfg.K_values(:)';  % row vector
N_values = cfg.N_values(:)';  % row vector

% ---- Baseline (per-device, K=1) -----------------------------------------
base_cfg       = cfg;
base_cfg.mode  = 'per_device';
base_cfg.K     = 1;
base           = calc_tdma_cspuf_throughput(base_cfg);
baseline_gp    = base.payload_goodput_bps;

if cfg.verbose
    fprintf('[Baseline] Per-device K=1 goodput: %.3f bps\n', baseline_gp);
end

% ---- Sweep and collect results ------------------------------------------
results = struct('mode',{},'N',{},'K',{},'goodput_bps',{},'improvement_pct',{});

% Per-device reuse
for K = K_values
    c        = cfg; 
    c.mode   = 'per_device';
    c.K      = K;
    met      = calc_tdma_cspuf_throughput(c);
    impr_pct = 100*(met.payload_goodput_bps / baseline_gp - 1);
    results(end+1) = struct('mode',"per_device",'N',NaN,'K',K, ...
        'goodput_bps',met.payload_goodput_bps,'improvement_pct',impr_pct); %#ok<AGROW>
end

% Broadcast once per frame (vary N and K)
for N = N_values
    for K = K_values
        c        = cfg;
        c.mode   = 'broadcast';
        c.K      = K;
        c.N      = N;
        met      = calc_tdma_cspuf_throughput(c);
        impr_pct = 100*(met.payload_goodput_bps / baseline_gp - 1);
        results(end+1) = struct('mode',"broadcast",'N',N,'K',K, ...
            'goodput_bps',met.payload_goodput_bps,'improvement_pct',impr_pct); %#ok<AGROW>
    end
end

% Convert to table
out.results = struct2table(results);

% ---- Plot 1: Per-device reuse -------------------------------------------
fig1 = figure('Name','Per-device challenge reuse');
pd = out.results(strcmp(out.results.mode,'per_device'),:);
plot(pd.K, pd.improvement_pct, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('Challenge reuse factor K (per-device)');
ylabel('Throughput improvement over baseline (%)');
title(sprintf('%s: Per-device challenge reuse', cfg.title_prefix));
xticks(unique(pd.K));
out.fig_per_device = fig1;

% ---- Plot 2: Broadcast once per frame -----------------------------------
fig2 = figure('Name','Broadcast once per frame');
hold on;
colsN = unique(out.results.N(~isnan(out.results.N)))';
for N = colsN
    row = out.results(strcmp(out.results.mode,'broadcast') & out.results.N==N, :);
    plot(row.K, row.improvement_pct, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'DisplayName', sprintf('N=%d devices', N));
end
hold off;
grid on;
xlabel('Challenge reuse factor K (broadcast per frame)');
ylabel('Throughput improvement over baseline (%)');
title(sprintf('%s: Broadcast challenge once per frame', cfg.title_prefix));
legend('Location','best');
xticks(unique(pd.K));
out.fig_broadcast = fig2;

% ---- Optional saves ------------------------------------------------------
if cfg.save_figs
    if ~exist(cfg.save_dir, 'dir'), mkdir(cfg.save_dir); end
    saveas(fig1, fullfile(cfg.save_dir, 'per_device_reuse.png'));
    saveas(fig2, fullfile(cfg.save_dir, 'broadcast_reuse.png'));
end

if cfg.save_csv
    if ~exist(cfg.save_dir, 'dir'), mkdir(cfg.save_dir); end
    writetable(out.results, fullfile(cfg.save_dir, 'lorawan_dr3_throughput_improvement.csv'));
end

% Return some extras
out.base_goodput_bps = baseline_gp;
out.cfg = cfg;

% -------------------- Local helpers --------------------------------------
function s = set_default(s, field, val)
    if ~isfield(s, field) || isempty(s.(field))
        s.(field) = val;
    end
end

end % <-- end main function

% ========== Helper: slot-time & goodput calculator =======================
function out = calc_tdma_cspuf_throughput(cfg)
% Compute TDMA payload goodput with explicit challenge overhead.
% Required: Rd, Ru, Lc, Lr, Lp, Tg, K, mode ('per_device'|'broadcast')
% Optional: N, Tpreamble_dl, Tpreamble_ul, hdr_dl_bits, hdr_ul_bits,
%           fec_dl_rate, fec_ul_rate

% Safe getter (prevents "Unrecognized field name" errors)
getf = @(S,f,d) (isfield(S,f) && ~isempty(S.(f))) * 1; %#ok<NASGU>
% (Use explicit if/else so MATLAB doesn't touch missing fields)
Rd           = ifget(cfg,'Rd',           1760.0);
Ru           = ifget(cfg,'Ru',           1760.0);
Lc           = ifget(cfg,'Lc',            256);
Lr           = ifget(cfg,'Lr',            128);
Lp           = ifget(cfg,'Lp',            512);
Tg           = ifget(cfg,'Tg',          5e-3);
Tpreamble_dl = ifget(cfg,'Tpreamble_dl',  0.0);
Tpreamble_ul = ifget(cfg,'Tpreamble_ul',  0.0);
hdr_dl_bits  = ifget(cfg,'hdr_dl_bits',   0);
hdr_ul_bits  = ifget(cfg,'hdr_ul_bits',   0);
fec_dl_rate  = ifget(cfg,'fec_dl_rate',   1.0);
fec_ul_rate  = ifget(cfg,'fec_ul_rate',   1.0);
K            = ifget(cfg,'K',             1);
mode         = lower(string(ifget(cfg,'mode','per_device')));
N            = ifget(cfg,'N',             1);

if K <= 0, error('K must be >= 1'); end
if any(strcmp(mode, "broadcast")) && (N < 1)
    error('Broadcast mode requires cfg.N >= 1');
end

% Bit->time (accounting for code rate)
bt_dl = @(bits) (bits ./ (Rd * fec_dl_rate));
bt_ul = @(bits) (bits ./ (Ru * fec_ul_rate));

% One DL challenge emission (full)
T_dl_once = bt_dl(Lc + hdr_dl_bits) + Tpreamble_dl;

% Amortized DL time per slot
switch mode
    case "per_device"
        T_dl_per_slot = T_dl_once / K;
    case "broadcast"
        T_dl_per_slot = T_dl_once / (N * K);
    otherwise
        error('Unknown mode: %s', mode);
end

% UL response and payload (always each slot)
T_ul_resp = bt_ul(Lr + hdr_ul_bits) + Tpreamble_ul;
T_ul_pay  = bt_ul(Lp);

% Total slot time and goodput
T_slot = T_dl_per_slot + T_ul_resp + T_ul_pay + Tg;

out = struct();
out.T_dl_per_slot        = T_dl_per_slot;
out.T_ul_response        = T_ul_resp;
out.T_ul_payload         = T_ul_pay;
out.T_guard              = Tg;
out.T_slot               = T_slot;
out.payload_goodput_bps  = Lp / T_slot;
out.fractions = struct( ...
   'challenge', T_dl_per_slot / T_slot, ...
   'response',  T_ul_resp      / T_slot, ...
   'payload',   T_ul_pay       / T_slot, ...
   'guard',     Tg             / T_slot);

% ---- nested safe getter ----
function v = ifget(S, f, d)
    if isfield(S, f) && ~isempty(S.(f))
        v = S.(f);
    else
        v = d;
    end
end

end
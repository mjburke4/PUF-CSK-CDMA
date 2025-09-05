function metrics = metrics_update(metrics, rx_bits, tx_bits, FAR, FRR, success, crps_used)
% metrics_update - Accumulate per-slot metrics into a struct
% metrics fields (per SNR index):
%   .bits_total, .bit_errors, .FAR_sum, .FRR_sum, .success_count, .slots, .crps_total

    if isempty(metrics)
        metrics.bits_total = 0;
        metrics.bit_errors = 0;
        metrics.FAR_sum = 0;
        metrics.FRR_sum = 0;
        metrics.success_count = 0;
        metrics.slots = 0;
        metrics.crps_total = 0;
    end

    metrics.bits_total    = metrics.bits_total + numel(tx_bits);
    metrics.bit_errors    = metrics.bit_errors + sum(rx_bits ~= tx_bits);
    metrics.FAR_sum       = metrics.FAR_sum + FAR;
    metrics.FRR_sum       = metrics.FRR_sum + FRR;
    metrics.success_count = metrics.success_count + double(success);
    metrics.slots         = metrics.slots + 1;
    metrics.crps_total    = metrics.crps_total + crps_used;
end
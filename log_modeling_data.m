function log_modeling_data(csv_path, scheme, node_id, snr_db, C_seed, obs)
% log_modeling_data - Append per-slot observations for offline modeling
%   csv_path: file path for CSV (created if missing)
%   scheme: 'CDMA' or 'CSK'
%   node_id: numeric
%   snr_db: scalar
%   C_seed: 1 x challenge_len (0/1) vector
%   obs: struct returned by txrx_*_node (contains rx waveform and stats)
%
% Writes one row per symbol (CSK) or per bit (CDMA).

    if ~exist('csv_path','var') || isempty(csv_path)
        csv_path = 'modeling_data.csv';
    end

    is_new = ~isfile(csv_path);
    fid = fopen(csv_path, 'a');

    if strcmpi(scheme, 'CDMA')
        if is_new
            fprintf(fid, 'scheme,node_id,snr_db,bit_idx,corr,conf,true_bit,rx_bit,seed_hash\n');
        end
        L = obs.L; Nbits = obs.Nbits;
        seed_hash = sum((1:numel(C_seed)) .* (double(C_seed(:))' + 1));
        for i = 1:Nbits
            true_bit = obs.tx_bits(i);
            rx_bit   = obs.corr_vals(i) > 0;
            fprintf(fid, 'CDMA,%d,%.3f,%d,%.6g,%.6g,%d,%d,%d\n', node_id, snr_db, i, obs.corr_vals(i), obs.conf_vals(i), true_bit, rx_bit, seed_hash);
        end
    else
        if is_new
            fprintf(fid, 'scheme,node_id,snr_db,sym_idx,maxc,second,correct,chosen,seed_hash\n');
        end
        Ns = obs.Ns;
        seed_hash = sum((1:numel(C_seed)) .* (double(C_seed(:))' + 1));
        for s = 1:Ns
            fprintf(fid, 'CSK,%d,%.3f,%d,%.6g,%.6g,%.6g,%d,%d\n', node_id, snr_db, s, obs.corr_max(s), obs.corr_2nd(s), obs.corr_true(s), obs.k_hat(s), seed_hash);
        end
    end
    fclose(fid);
end
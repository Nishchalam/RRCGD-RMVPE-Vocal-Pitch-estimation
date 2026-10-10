%% Run Root Extraction — PTDB-TUG, 32ms window, 10ms hop, 8kHz
% Run from inside this directory after setup_lf_path / rootsIni (see
% this directory's README.md). Edit DATA_ROOT below to point at the
% output of data_pipeline/3_augment/generate_noisy_data.py (i.e. the
% dataset/PTDB_data_10ms_hop/ layout produced by steps 1-3).

DATA_ROOT = '../../dataset/PTDB_data_10ms_hop';

for split = ["train", "valid", "test"]
    data_path = fullfile(DATA_ROOT, split);
    fprintf('PTDB 32ms window - %s\n', split);
    extract_roots_batch(data_path, 8000, 32, 10);
end

% To extract a different window length (e.g. for the 64ms-window
% variant some configs support), just change the 32 above -- the output
% subdirectory name (roots_WL_<ms>ms_HOP_<ms>ms) is derived automatically.

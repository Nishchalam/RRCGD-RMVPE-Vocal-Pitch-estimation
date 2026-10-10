# Data pipeline: raw audio → trained model → results

End-to-end steps to go from nothing to reproduced results, each in its
own numbered subfolder.

| Step | What | Where |
|---|---|---|
| 1 | Download PTDB-TUG + set up the Python environment | [1_download/](1_download/) |
| 2 | Build the speaker-disjoint train/valid/test split | [2_split/](2_split/) |
| 3 | Generate SNR-augmented noisy data (for clean+noisy training) | [3_augment/](3_augment/) |
| 4 | Extract polynomial roots (MATLAB) | [4_roots_matlab/](4_roots_matlab/) |
| 5 | Edit a YAML config for your experiment list | see below |
| 6 | Run `run_experiments.py` | see below |
| 7 | Aggregate and analyze results | see below |

Run them in order — each step's output is the next step's input.

```bash
# 1. Download + env setup
cd data_pipeline/1_download
./download_ptdb.sh /path/to/PTDB_raw
conda env create -f ../../environment.yml -n gd_exps && conda activate gd_exps
pip install torch==2.10.0 torchaudio==2.10.0 torchvision==0.25.0 --index-url https://download.pytorch.org/whl/cu128
pip install -r ../../requirements.txt

# 2. Speaker-disjoint split
cd ../2_split
python split_ptdb.py --raw-dir /path/to/PTDB_raw --out ../../dataset/PTDB_data_10ms_hop

# 3. Noisy augmentation (skip if you only want clean-only configs)
cd ../3_augment
python generate_noisy_data.py --dataset-root ../../dataset/PTDB_data_10ms_hop

# 4. Root extraction (MATLAB -- see 4_roots_matlab/README.md for full instructions)
cd ../4_roots_matlab
# in MATLAB: setup_lf_path; rootsIni; run_experiment;

# 5. Pick/edit a config
cd ../..
#   configs/experiments_clean_256bins_bs8.yaml         (clean-only, 256 bins)
#   configs/experiments_clean_512bins_bs8.yaml         (clean-only, 512 bins)
#   configs/experiments_clean_noisy_256bins_bs8.yaml   (clean+noisy, 256 bins)
#   configs/experiments_clean_noisy_512bins_bs8.yaml   (clean+noisy, 512 bins)
# each sweeps 7 techniques x 23 alphas -- see the main README for what
# E1-E4/mel/linear mean, or --only one experiment name to run a subset

# 6. Dry run, then train
python run_experiments.py --config configs/experiments_clean_256bins_bs8.yaml --dry-run
python run_experiments.py --config configs/experiments_clean_256bins_bs8.yaml --skip-existing

# 7. Aggregate results
python export_validation_csv.py
python compile_best_results.py
python compile_results_by_condition.py   # per-SNR breakdown, clean+noisy runs only
```

Full details, caveats, and reproducibility notes for the training/eval
side (step 5 onward) are in the top-level [INSTALL.md](../INSTALL.md)
and [README.md](../README.md).

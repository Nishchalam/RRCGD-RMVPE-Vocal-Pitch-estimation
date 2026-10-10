# Step 3: SNR-augmented noisy data (for clean+noisy training)

```bash
python generate_noisy_data.py --dataset-root dataset/PTDB_data_10ms_hop
```

Adds white Gaussian noise at 4 target SNR levels (20, 10, 5, 0 dB) to
every clean `.wav` in `train/`, `valid/`, and `test/` from [step
2](../2_split/), writing `snr{20,10,05,00}_<stem>.wav` next to the
original. Seeded (`RANDOM_SEED = 42`), so this is reproducible. A
`noise_generation_log.csv` records the target vs. actually-measured
SNR per generated file.

This step is **required** for the `configs/experiments_clean_noisy_*.yaml`
configs (clean+noisy training, this project's best-performing setting)
and **optional** if you only want to run the clean-only configs
(`configs/experiments_clean_*.yaml` without `_noisy`).

No label files are generated for the noisy wavs — `src/dataset.py`
resolves `snrNN_<stem>`'s pitch label back to the clean `<stem>.npy`
by stripping the SNR prefix, so there's nothing to duplicate.

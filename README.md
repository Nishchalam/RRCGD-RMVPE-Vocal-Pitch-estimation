# RRCGD-RMVPE: Vocal Pitch Estimation

Pitch (F0) estimation combining **RRCGD** (Reflected Root Chirp Group
Delay spectrum), a feature from Singh, Batra & Ramalingam (IIT Madras,
ICASSP 2025), with **RMVPE**, an established pitch-estimation
architecture (DeepUnet → BiGRU → salience head), in place of RMVPE's
native mel-spectrogram input.

This repository contains the training/evaluation code and
representative experiment configs. It does **not** contain any audio
data, labels, or trained checkpoints — see [INSTALL.md](INSTALL.md) for
how to set up the environment and point it at your own copy of
PTDB-TUG.

## Summary

On PTDB-TUG (speech), RRCGD→RMVPE reaches **94.75% RPA** (50-cent
threshold, leakage-free speaker-disjoint split, clean+noisy training —
see below), exceeding every number in the RRCGD paper's own Table II —
ACS-of-RRCGD 92.51%, plain RRCGD 91.91%, CREPE 88.55% — obtained under
that paper's original Conv-LSTM + Transformer architecture. This
isolates one question: swapping in RMVPE's architecture, holding the
RRCGD-family feature fixed, improves pitch estimation on speech.

## Results

**Best result per completed sweep** (bs=8, the adopted standard; α =
the focal-loss weighting term; all rows 256 bins, 23/23 alphas unless
noted):

| Training data | Best technique | Test RPA | α |
|---|---|---|---|
| **Clean+noisy, 256 bins** | **E3-linear** (upper-freq-truncated) | **94.75%** | 18.5 |
| Clean+noisy, 256 bins | E4-mel (index-truncated) | 94.73% | 19.0 |
| Clean+noisy, 512 bins | E3-linear (upper-freq-truncated) | 94.66% | 20.0 |
| Clean-only, 256 bins | E4-mel (index-truncated) | 94.51% | 20.0 |
| Clean-only, 512 bins | E4-mel | 93.92% | 16.5 |

**Clean-only vs. clean+noisy training, evaluated on identical noisy
test data** (mean over 7 techniques, 256 bins, bs=8):

| Test condition | Clean-only training | Clean+noisy training | Δ |
|---|---|---|---|
| Clean | 94.18% | 94.56% | +0.37 |
| 20dB SNR | 94.11% | 94.61% | +0.50 |
| 10dB SNR | 93.29% | 94.56% | +1.27 |
| 5dB SNR | 91.57% | 94.24% | +2.67 |
| 0dB SNR | 86.96% | 93.38% | **+6.43** |
| Overall (frame-weighted) | 92.02% | 94.27% | +2.25 |

A clean-only model loses ~7.2 points of RPA from its own clean-condition
score at 0dB SNR. Training on clean + all four SNR levels together
closes nearly all of that gap and, incidentally, improves the
clean-condition score too — consistent with the added data acting as a
regularizer rather than trading clean accuracy for robustness. The best
clean+noisy technique (94.75%) exceeds the clean-only headline even on
the clean-condition metric alone, not just on robustness.

## Method

**Data.** PTDB-TUG, 8kHz, 10ms hop, speaker-disjoint train/test split.

**Feature.** RRCGD replaces each frame's mel-spectrogram with a
group-delay spectrum computed from the reflected roots of the frame's
samples treated as a polynomial (`numpy.roots`), evaluated on a
frequency grid. Four grid variants are swept per bin family (mel,
linear):

- **E1**: full range to `fs/2`. **E2**: full range to `fs/4`.
- **E3**: build the full `0→fs/2` grid, then keep only bins with center
  frequency `≤ fs/4` ("upper-frequency-truncated").
- **E4**: build the full grid, then keep the first half of bins **by
  index**, not frequency ("upper-bins-truncated") — on a linear grid
  this is identical to E3, so only the mel variant is run.

Implementation: [`src/spec.py`](src/spec.py), `RRCGDSpectrogram`.

**Training.** α (focal-loss weighting) is swept over 9.0–20.0 per
technique; batch size and bin count (256 vs. 512) are swept
independently. Training runs in fp32 — fp16 autocast caused some
high-α runs (mostly 512-bin, full-range techniques) to collapse to a
non-finite loss with no checkpoint; fp32 does not. The ±50 clamp inside
the RRCGD sum (`src/spec.py`) is load-bearing for the same reason — see
[INSTALL.md](INSTALL.md) §7.

## Findings

**Batch size.** bs=8 has the highest median RPA across the alpha sweep
at both 256 and 512 bins (not merely the best single alpha); bs=32
collapses to 87–90% RPA on a subset of alphas, most severely at 512
bins.

**Bin count.** 256 bins edges out 512 bins technique-for-technique at
matched batch size; the margin is consistently small.

**Noise robustness.** See Results above.

## Getting started

```bash
# 1. Set up the environment and point it at your own PTDB-TUG copy
#    (see INSTALL.md for full details)

# 2. See every (experiment, alpha) job a config would run, without running it
python3 run_experiments.py --config configs/experiments_clean_256bins_bs8.yaml --dry-run

# 3. Run one experiment's full alpha grid; resumable
python3 run_experiments.py --config configs/experiments_clean_256bins_bs8.yaml \
  --only ptdb_10ms_clean_linear_fs2 --skip-existing

# 4. Aggregate results across the sweep
python3 export_validation_csv.py
python3 compile_best_results.py
```

Full setup, data layout, and reproducibility notes:
[INSTALL.md](INSTALL.md). For the complete pipeline starting from
nothing — downloading PTDB-TUG, building the speaker-disjoint split,
generating noisy data, and extracting roots in MATLAB — see
[data_pipeline/](data_pipeline/).

## Repository layout

```
data_pipeline/                # steps 1-4: raw audio -> trained-model-ready data (see data_pipeline/README.md)
  1_download/                   # PTDB-TUG download + environment setup
  2_split/                       # speaker-disjoint train/valid/test split
  3_augment/                     # SNR-augmented noisy data (for clean+noisy training)
  4_roots_matlab/                 # polynomial root extraction (MATLAB, Lindsey-Fox toolbox)
run_experiments.py           # orchestrates one run per (experiment, alpha) in a config
train_single.py               # single training run; reads EXP_* env vars run_experiments.py sets
evaluate.py                   # generic evaluation entry point
evaluate_ptdb_checkpoint.py    # scores a checkpoint against the PTDB test split (per-SNR breakdown)
export_validation_csv.py      # aggregates every run dir's results into one CSV
compile_best_results.py       # picks the best alpha per technique from that CSV
compile_results_by_condition.py   # per-SNR-condition RPA breakdown for clean+noisy runs
check_repro_env.py            # environment/determinism fingerprint for cross-machine debugging
src/
  model.py, deepunet.py        # RMVPE architecture
  spec.py                      # RRCGDSpectrogram -- the RRCGD feature extractor
  dataset.py                    # PTDB dataset loading, SNR-condition handling, frame caching
  loss.py, inference.py, seq.py, utils.py, constants.py
configs/
  experiments_clean_256bins_bs8.yaml         # clean-only, 256 bins
  experiments_clean_512bins_bs8.yaml         # clean-only, 512 bins
  experiments_clean_noisy_256bins_bs8.yaml   # clean+noisy, 256 bins
  experiments_clean_noisy_512bins_bs8.yaml   # clean+noisy, 512 bins
```

## Reference

Singh, N., Batra, M., & Ramalingam, C. S. "Improved Pitch and Voicing
Determination Using the Reflected Root Chirp Group Delay Spectrum."
ICASSP 2025. DOI: 10.1109/ICASSP49660.2025.10890591.

## License

No license file is included yet — all rights reserved by default until
one is added. Add a `LICENSE` file before relying on this repo being
reusable by others.

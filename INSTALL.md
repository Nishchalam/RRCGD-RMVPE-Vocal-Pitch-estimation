# Installation & Reproducibility Guide

How to set up this project on a new machine, and exactly what
reproducibility you should and shouldn't expect once you do.

## 1. What "reproducible" means here — read this first

This codebase seeds every RNG it touches (`src/utils.py`'s `set_seed()`)
and forces PyTorch/cuDNN/cuBLAS into deterministic-kernel mode when
`deterministic: true` (the project default). That gets you:

- **Bit-exact reproduction**, if you match: same GPU **architecture**
  (not just vendor — an RTX 4090 and an A100 are different), same
  PyTorch/CUDA/cuDNN/driver versions, same OS. This is the strongest
  guarantee PyTorch can make today — nobody's determinism setup does
  better than "same hardware + same software stack."
- **Statistically equivalent, not bit-identical**, results across
  *different* GPU architectures, even with every version pinned
  identically. Floating-point addition isn't associative, and different
  GPU architectures reduce sums (matmuls, convolutions, this project's
  own chunked RRCGD accumulation) in a different order — the gap is a
  few ULP per op but compounds over a full training run into a
  real-but-small difference in the final RPA. Expect numbers within
  noise of each other, not identical to the last decimal.

So: **pin your GPU model and software stack (this guide) if you need
bit-exact reproduction.** If you're on different hardware, match the
software stack anyway and treat "reproducible" as "same result within
normal run-to-run noise," not "bit-identical."

## 2. Reference environment (what the original results were trained on)

| | |
|---|---|
| GPU | NVIDIA GeForce RTX 4090 (compute capability 8.9) |
| OS | Ubuntu 22.04 LTS |
| Python | 3.10.18 |
| PyTorch | 2.10.0+cu128 |
| CUDA (via PyTorch) | 12.8 |
| cuDNN (via PyTorch) | 9.10.2 |

`requirements.txt` is the curated list of what this codebase actually
imports.

## 3. Set up the environment

```bash
# 1. Create the Python environment
conda env create -f environment.yml -n rrcgd_rmvpe
conda activate rrcgd_rmvpe

# 2. Install torch/torchaudio/torchvision from PyTorch's own CUDA wheel
#    index -- NOT plain PyPI/conda-forge, which can silently give you a
#    CPU-only or mismatched-CUDA build with the same version string.
pip install torch==2.10.0 torchaudio==2.10.0 torchvision==0.25.0 \
    --index-url https://download.pytorch.org/whl/cu128

# 3. Everything else
pip install -r requirements.txt
```

**Different CUDA version on your machine?** Pick the matching PyTorch
wheel index instead of cu128 (e.g. `cu121`, `cu124`) — check
https://pytorch.org/get-started/locally/ for the current index URLs.
This changes the bit-exactness guarantee (§1): a different CUDA/cuDNN
build is not bit-identical to the reference runs, even on identical
hardware.

### Verify the install

```bash
python3 -c "
import torch
print('torch', torch.__version__)
print('cuda available:', torch.cuda.is_available())
print('cuda version:', torch.version.cuda)
print('cudnn version:', torch.backends.cudnn.version())
print('gpu:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'none')
print('compute capability:', torch.cuda.get_device_capability(0) if torch.cuda.is_available() else 'n/a')
"
```

Compare the output against §2. Everything matching except the GPU name
is fine as long as compute capability matches too, if you're after
bit-exactness; a different compute capability means treat results as
statistically-equivalent, not bit-identical.

`check_repro_env.py` writes a fuller `env_report_<hostname>.json`
fingerprint if you want to diff two machines directly:

```bash
python3 check_repro_env.py
```

## 4. Get the dataset

This repo contains no audio data or labels — only code. You need
**PTDB-TUG** (speech, 8kHz) yourself:

- PTDB-TUG is a published academic speech corpus (Pitch Tracking
  Database from Graz University of Technology); obtain it from its
  official distribution.
- This codebase expects the data **already split into speaker-disjoint
  `train/`/`test/` folders** with each utterance's pitch reference
  already extracted, laid out as:

  ```
  dataset/
    PTDB_data_10ms_hop/
      train/  test/
        <stem>.wav   <stem>.f0   <stem>.npy    # audio + pitch label
        roots_WL_32ms_HOP_10ms/                # pre-extracted complex roots, 32ms window @ 8kHz
        _frame_cache_WL*ms_HOP10ms_FREF10.0_NCLASS*/   # auto-generated on first run, safe to delete
  ```

- The `roots_WL_<window_samples>ms_HOP_<hop>ms/` directories hold
  pre-extracted complex polynomial roots per frame (see
  `src/spec.py`'s `RRCGDSpectrogram` for how they're consumed) — this
  is a one-time preprocessing step over the raw audio, not something
  `train_single.py` does itself. You'll need to write that extraction
  step against your own copy of PTDB-TUG (root-finding via
  `numpy.roots` on each windowed frame) before training.
- `_frame_cache_*` directories are a pure feature cache
  (`src/dataset.py`'s `_ensure_cache`), auto-created on first access —
  safe to delete anytime to reclaim disk; they just regenerate.

Point `data_paths.ptdb` in whichever `configs/experiments_*.yaml` you
use at your local copy of this layout.

## 5. Running an experiment

```bash
# Dry run first -- prints every job it would launch without running anything
python3 run_experiments.py --config configs/experiments_clean_256bins_bs8.yaml --dry-run

# One experiment, full alpha grid (resumable: a killed/crashed run just
# picks back up from its last completed alpha next time)
python3 run_experiments.py --config configs/experiments_clean_256bins_bs8.yaml \
  --only ptdb_10ms_clean_linear_fs2 --skip-existing

# Reproducibility-relevant CLI overrides (see run_experiments.py --help
# for the full list): --seed, --deterministic {true,false},
# --batch-size, --patience-epochs, --use-cache {true,false}, --no-seed
```

The four provided configs cover the project's main grid:

| Config | Bins | Training data |
|---|---|---|
| `experiments_clean_256bins_bs8.yaml` | 256 | clean only |
| `experiments_clean_512bins_bs8.yaml` | 512 | clean only |
| `experiments_clean_noisy_256bins_bs8.yaml` | 256 | clean + 4 SNR levels |
| `experiments_clean_noisy_512bins_bs8.yaml` | 512 | clean + 4 SNR levels |

Each sweeps 7 techniques (E1–E4 frequency representation × {mel,
linear} bin spacing, with E4 mel-only) × 23 alpha values — see the main
README for what E1–E4 mean.

## 6. Aggregating results

```bash
python3 export_validation_csv.py     # one row per (experiment, alpha) checkpoint
python3 compile_best_results.py      # best alpha per technique
python3 compile_results_by_condition.py   # per-SNR-condition RPA breakdown (clean+noisy runs only)
```

## 7. Known gotchas (learned the hard way)

- **GPU memory headroom from a single `nvidia-smi` snapshot is
  unreliable** when multiple jobs share a card — a snapshot showing
  free memory can still OOM the next launch. Budget per-job memory
  empirically rather than trusting a live snapshot, and leave real
  margin.
- **Batch size preference is dataset- and technique-specific**, not a
  fixed architectural property — verified independently at multiple
  bin counts for PTDB (bs=8 wins here) — check directly for any new
  dataset/config rather than assuming this transfers.
- **`src/spec.py`'s `_chunked_rrcgd_sum` clamp is load-bearing.** It
  keeps a numerically-unstable term (`num/den` can spike when a root is
  near-resonant with a grid frequency) inside a range fp16 autocast can
  handle. Removing it causes real non-finite-loss training instability
  within a few dozen alphas — don't remove it.

import os

# ── Sample rate ───────────────────────────────────────────────────────────────
# MIR-1K (music) → 16 000 Hz  (default)
# PTDB   (speech) → 8 000 Hz  — set EXP_SAMPLE_RATE=8000 via run_experiments.py
SAMPLE_RATE = int(os.environ.get("EXP_SAMPLE_RATE", 16000))

N_CLASS  = 360
N_BINS   = int(os.environ.get("EXP_N_BINS",256))
R0       = 1.008

# ── Overridable per-experiment via environment variables ──────────────────────
# run_experiments.py sets EXP_* before spawning each subprocess so every fresh
# Python process sees the correct values before any `from src import …` fires.
# Fallbacks reproduce the original hard-coded behaviour.
F_MAX         = int(os.environ.get("EXP_F_MAX",         SAMPLE_RATE // 2))
WINDOW_LENGTH = int(os.environ.get("EXP_WINDOW_LENGTH", 2048))
HOP_MS        = int(os.environ.get("EXP_HOP_MS", 20))

# Base bin methods:
# "mel_aligned" → Mel-frequency-aligned grid from 0 to F_MAX
# "linear_fmax" → linear-frequency grid from 0 to F_MAX
#
# Truncation variants:
# "*_upper_frequency_truncated"
#     Build the original grid from 0 to fs/2, then retain bins <= F_MAX.
#
# "*_upper_bins_truncated"
#     Build the original grid from 0 to fs/2, then retain the lower
#     half of the bins. F_MAX identifies the experiment cutoff but
#     does not determine bin selection.
BIN_METHOD = os.environ.get("EXP_BIN_METHOD", "mel_aligned")

# PTDB data selection / caching. These affect only dataset I/O, not the model.
CLEAN_ONLY = os.environ.get("EXP_CLEAN_ONLY", "0") == "1"
USE_CACHE = os.environ.get("EXP_USE_CACHE", "1") == "1"
N_DENSE    = int(os.environ.get("EXP_N_DENSE", "2048"))  # dense grid size before pooling

CONST    = 1997.3794084376191
PITCH_TH = float(os.environ.get("EXP_PITCH_TH", 0.5))

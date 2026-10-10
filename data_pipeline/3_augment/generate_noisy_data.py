"""
Step 3: Generate SNR-augmented noisy PTDB data (clean + noisy training).

Input structure (the output of step 2 -- note: labels are .npy here,
not .f0, since step 2 already converted them):
PTDB_data_10ms_hop/
    train/
        xxx.wav
        xxx.npy
    valid/
        ...
    test/
        ...

Output
------
The original files remain untouched. Only noisy .wav files are
generated -- NOT a copy of the label. src/dataset.py's PTDB class
resolves a noisy root/checkpoint back to its clean .npy label by
stripping the SNR prefix (strip_snr_prefix()), so a label never needs
to exist under the snrNN_ prefix at all.

Generated files:
    snr20_xxx.wav
    snr10_xxx.wav
    snr05_xxx.wav
    snr00_xxx.wav

Usage:
    python generate_noisy_data.py --dataset-root dataset/PTDB_data_10ms_hop
"""

import argparse
from pathlib import Path
import numpy as np
import soundfile as sf
import pandas as pd
from tqdm import tqdm

# ==========================================================
# Configuration
# ==========================================================

SPLITS = [
    "train",
    "valid",
    "test",
]

SNR_LEVELS = [
    20,
    10,
    5,
    0,
]

RANDOM_SEED = 42

rng = np.random.default_rng(RANDOM_SEED)

# ==========================================================
# Functions
# ==========================================================


def signal_power(x):
    return np.mean(x.astype(np.float64) ** 2)


def add_awgn(signal, snr_db):
    """
    Add white Gaussian noise to achieve the desired SNR.
    """

    p_signal = signal_power(signal)

    p_noise = p_signal / (10 ** (snr_db / 10))

    noise = rng.standard_normal(signal.shape)

    noise = noise / np.sqrt(signal_power(noise))

    noise *= np.sqrt(p_noise)

    noisy = signal + noise

    return noisy.astype(signal.dtype)


def measured_snr(clean, noisy):
    noise = noisy - clean

    ps = signal_power(clean)
    pn = signal_power(noise)

    return 10 * np.log10(ps / pn)


# ==========================================================
# Main
# ==========================================================


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset-root", required=True, type=Path,
                     help="e.g. dataset/PTDB_data_10ms_hop (step 2's output)")
    args = ap.parse_args()
    dataset_root = args.dataset_root
    log_file = dataset_root / "noise_generation_log.csv"

    rows = []

    for split in SPLITS:

        split_dir = dataset_root / split

        # only augment original (non-snr-prefixed) files, in case this is re-run
        wav_files = sorted(p for p in split_dir.glob("*.wav") if not p.stem.startswith("snr"))

        print(f"\n{split}: {len(wav_files)} wav files")

        for wav_path in tqdm(wav_files):

            stem = wav_path.stem

            npy_path = wav_path.with_suffix(".npy")

            if not npy_path.exists():
                print(f"WARNING: Missing label {npy_path}")
                continue

            audio, sr = sf.read(wav_path)

            # Mono only
            if audio.ndim > 1:
                audio = audio.mean(axis=1)

            duration = len(audio) / sr

            for snr in SNR_LEVELS:

                noisy = add_awgn(audio, snr)

                out_name = f"snr{snr:02d}_{stem}"
                out_wav = split_dir / f"{out_name}.wav"

                sf.write(out_wav, noisy, sr)

                # No label is written for snrNN_<stem> -- src/dataset.py
                # resolves it back to <stem>.npy by stripping the prefix.

                actual_snr = measured_snr(audio, noisy)

                rows.append(
                    {
                        "split": split,
                        "original": wav_path.name,
                        "generated": out_wav.name,
                        "target_snr_db": snr,
                        "measured_snr_db": actual_snr,
                        "sample_rate": sr,
                        "duration_sec": duration,
                    }
                )

    # ==========================================================
    # Save log
    # ==========================================================

    df = pd.DataFrame(rows)
    df.to_csv(log_file, index=False)

    print("\nDone.")
    print(f"Generated {len(df)} noisy files.")
    print(f"Log written to:\n{log_file}")


if __name__ == "__main__":
    main()
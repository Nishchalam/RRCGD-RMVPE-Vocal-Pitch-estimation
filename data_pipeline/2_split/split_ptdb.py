"""
Step 2: Build a speaker-disjoint train/valid/test split of PTDB-TUG.

Input: the raw, unzipped PTDB-TUG download (step 1) -- a directory tree
containing, somewhere under it, paired files named
`mic_<speaker>_<utterance>.wav` (audio) and
`ref_<speaker>_<utterance>.f0` (pitch reference: whitespace-separated
columns, column 0 = pitch in Hz, 0 = unvoiced, one row per 10ms frame).
The exact subfolder nesting is not assumed -- this script searches
recursively.

Output: a flat layout this project's src/dataset.py expects:

    <out>/
      train/  <speaker>_<utt>.wav   <speaker>_<utt>.npy
      valid/  (same)
      test/   (same)

`.npy` is just the reference `.f0` file's pitch column (Hz), saved as a
1D float array -- src/dataset.py loads it directly as the per-frame
pitch target, frame-index-aligned with the roots that step 4 computes
from the matching `.wav` (both run at the project's native 10ms hop,
so no resampling/interpolation of the label is needed here).

This is this project's actual speaker split (12 train / 4 valid / 4
test, PTDB-TUG's 20 speakers): change SPEAKER_SPLIT below if you want
a different split, but keep it speaker-disjoint -- the whole point of
this project's results is that no speaker appears in more than one of
train/valid/test.

Usage:
    python split_ptdb.py --raw-dir /path/to/PTDB_raw --out dataset/PTDB_data_10ms_hop
    python split_ptdb.py --raw-dir /path/to/PTDB_raw --out dataset/PTDB_data_10ms_hop --dry-run
"""
import argparse
import shutil
from collections import defaultdict
from pathlib import Path

import numpy as np

SPEAKER_SPLIT = {
    "train": {"F01", "F02", "F03", "F04", "F05", "F06",
              "M01", "M02", "M03", "M04", "M05", "M06"},
    "valid": {"F07", "F08", "M07", "M08"},
    "test":  {"F09", "F10", "M09", "M10"},
}

SPEAKER_TO_SPLIT = {spk: split for split, spks in SPEAKER_SPLIT.items() for spk in spks}


def speaker_of(stem_after_prefix: str) -> str:
    """'F01_sa1' -> 'F01'"""
    return stem_after_prefix.split("_")[0]


def find_pairs(raw_dir: Path):
    """Recursively find every (wav, f0) pair under raw_dir, matched by
    stem after stripping the 'mic_'/'ref_' prefix PTDB-TUG uses."""
    wavs = {p.stem[len("mic_"):]: p for p in raw_dir.rglob("mic_*.wav")}
    refs = {p.stem[len("ref_"):]: p for p in raw_dir.rglob("ref_*.f0")}
    common = sorted(set(wavs) & set(refs))
    missing_ref = sorted(set(wavs) - set(refs))
    missing_wav = sorted(set(refs) - set(wavs))
    return [(wavs[s], refs[s], s) for s in common], missing_ref, missing_wav


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw-dir", required=True, type=Path,
                     help="Root of the extracted PTDB-TUG download (step 1's output)")
    ap.add_argument("--out", required=True, type=Path,
                     help="Output root, e.g. dataset/PTDB_data_10ms_hop")
    ap.add_argument("--dry-run", action="store_true",
                     help="Print what would be written without writing anything")
    args = ap.parse_args()

    pairs, missing_ref, missing_wav = find_pairs(args.raw_dir)
    print(f"Found {len(pairs)} matched (wav, f0) pairs under {args.raw_dir}")
    if missing_ref:
        print(f"  WARNING: {len(missing_ref)} wav file(s) with no matching .f0 "
              f"(e.g. {missing_ref[:3]})")
    if missing_wav:
        print(f"  WARNING: {len(missing_wav)} .f0 file(s) with no matching wav "
              f"(e.g. {missing_wav[:3]})")

    counts = defaultdict(int)
    skipped_speakers = set()

    for wav_path, f0_path, stem in pairs:
        speaker = speaker_of(stem)
        split = SPEAKER_TO_SPLIT.get(speaker)
        if split is None:
            skipped_speakers.add(speaker)
            continue

        out_dir = args.out / split
        out_wav = out_dir / f"{stem}.wav"
        out_npy = out_dir / f"{stem}.npy"

        if args.dry_run:
            print(f"  [DRY] {wav_path.name} -> {out_wav}")
            counts[split] += 1
            continue

        out_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(wav_path, out_wav)

        pitch = np.loadtxt(f0_path, usecols=0, ndmin=1).astype(np.float64)
        np.save(out_npy, pitch)

        counts[split] += 1

    print("\nSummary:")
    for split in ("train", "valid", "test"):
        print(f"  {split:>5}: {counts[split]} file(s)  "
              f"(speakers: {sorted(SPEAKER_SPLIT[split])})")
    if skipped_speakers:
        print(f"\n  Speakers found but not in SPEAKER_SPLIT (skipped): "
              f"{sorted(skipped_speakers)}")

    if args.dry_run:
        print("\nDry run -- nothing was written. Re-run without --dry-run to write files.")


if __name__ == "__main__":
    main()

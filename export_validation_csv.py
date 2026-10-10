#!/usr/bin/env python3
"""
Flatten all completed validation runs (result.txt) under runs_seeded/ into
one CSV -- the running "all results so far" export, matching the column
format used by validation_result.ipynb and the older gen_results_csv.py.

This is the PREDECESSOR step to compile_best_results.py: this dumps every
(experiment, alpha) row with its last validation checkpoint; that script
later reduces this down to one best-alpha row per experiment using
test-set inference. Safe to re-run any time -- it always reflects
whatever has completed so far under runs_seeded/.

Fired automatically every 10 completed jobs by auto_export_every_10.sh.

Sweep-constant fields (hop_ms, dataset, SR, data_condition, window_ms,
batch_size, N_Bins) don't vary per row in a single sweep, so they're
pulled out into the filename instead of repeating them on every row.
(A companion <csv-name>_meta.json used to carry these plus a few other
constants too, but was dropped 2026-09-19 -- once wl/bs/N_Bins became
explicit in every filename, meta.json was redundant with it.)

Filename shape: <dataset>_<hop>ms_sr<SR>_<condition>_wl<wl>ms_bs<bs>_<n>bins_validation_results.csv
-- always includes window-length/batch-size/bin-count explicitly (not just
when they differ from some "default"), since these CSVs are meant to be
directly plottable/publishable and a uniform naming scheme matters more
than a shorter one. See compile_best_results.py's matching output-filename
convention (best_results_<dataset>_<condition>_wl<wl>ms_bs<bs>_<n>bins.csv).
"""

import re
import csv
from datetime import datetime
from pathlib import Path

RUNS_DIR = Path(__file__).parent / "runs_seeded" / "pitch_th_0.5"
OUT_DIR = Path(__file__).parent

_SR_FL_BINS_NEW = re.compile(r"(\d+)_FL_(\d+)_bins(\d+)")
_SR_FL_BINS_OLD = re.compile(r"(\d+)_FL_(\d+)_(\d+)_bins")
# The trailing (?:_(seed\d+|noseed))? is 2026-09-22's fix: run_experiments.py's
# logdir_for() appends "_seed<N>" for a non-default seed or "_noseed" for
# --no-seed (see its own docstring), but this regex's .fullmatch() had no
# way to match either suffix -- any seed-tagged run directory silently hit
# the `if not m: skipped += 1; continue` branch below and vanished from
# every CSV, with zero indication anything was dropped. This is exactly
# what happened to 2026-09-19's seed=123 ablation: its rows only ever
# survived in raw *_cleanonly_logs/*.log files, never any committed CSV.
# No seed-tagged directories exist on disk as of this fix (nothing to
# recover), but it's latent-bug prevention for the next one.
_FEATURE_DIR = re.compile(r"(.+)_([0-9]*\.?[0-9]+)_(\d+)(?:_(seed\d+|noseed))?")
_PITCH_DIR = re.compile(r"Pitch_FL_([0-9]*\.?[0-9]+)_([0-9]*\.?[0-9]+)")
# Matches the dataset_tag segment train_single.py writes, e.g. "ptdb_sr8000"
# or "mir1k_sr16000" -- see train_single.py's `dataset_tag = f"{DATASET}_sr{SAMPLE_RATE}"`.
_DATASET_SR = re.compile(r"^([a-zA-Z0-9]+)_sr(\d+)$")

# hop_ms isn't recorded anywhere in the run directory path itself (unlike
# N_Bins/F_MAX/window_length), so it has to come from a lookup -- these are
# each dataset's current fixed sweep value (see the datasets' respective
# experiments_*.yaml `defaults.hop_length`).
HOP_MS_BY_DATASET = {"ptdb": 10, "mir1k": 20}

# window_length (FL, in samples) maps to a window duration in ms that's
# constant across a whole sweep but NOT across sweeps -- MIR-1K now has both
# a 512-sample/32ms track and a 2048-sample/128ms track running concurrently,
# and PTDB now has multiple batch sizes running concurrently at the same
# window length (256bs8/256bs16/256bs32/512bs8/512bs16/512bs32, 2026-09-19 --
# caught the 256-bin sweep silently running at bs=32 instead of the
# project's finalized bs=8). Without folding window_ms AND batch_size into
# the grouping key below, rows from incompatible training configurations
# would land in the same *_validation_results.csv (sharing only dataset +
# data_condition), and compile_best_results.py would pick a "best alpha" by
# RPA across them as if they were one experiment.
#
# 2026-09-19: switched from "tag only if non-default" to ALWAYS including
# wl/bs/n_bins explicitly in every filename (user request -- these CSVs
# feed plotting/publishing directly, and a uniform, fully-self-describing
# naming scheme beats one where half the files carry no tag because they
# happen to match today's default). This renames every existing PTDB CSV;
# see this directory's README.md for the current filenames.

# Per-row columns only -- sweep-constant fields live in the meta JSON instead.
# RPA/RCA/OA/VR/VFA = validation (best checkpoint, from result.txt).
# test_RPA/etc = held-out test-set inference on that same best checkpoint
# (test_eval.csv, "clean" condition) -- blank until eval_test_set() has
# run for that alpha. This is what compile_best_results.py picks the best
# alpha by.
# Also includes separate columns for each noise condition (snr00, snr05, snr10, snr20)
# so we can see how the model performs across the full noise spectrum.
COLUMNS = [
    "FL", "N_Bins", "F_MAX", "feature_type", "technique",
    "batch_size", "seed_tag",
    "alpha", "gamma", "iteration",
    "RPA", "RCA", "OA", "VR", "VFA",
    "test_RPA", "test_RCA", "test_OA", "test_VR", "test_VFA",
    "test_RPA_snr00", "test_RPA_snr05", "test_RPA_snr10", "test_RPA_snr20",
]

# (dataset, feature_type, F_MAX) -> E1-E4, per each dataset's experiments_*.yaml
# representations. F_MAX values collide across datasets at their native rate
# (MIR-1K's E2 F_MAX=4000 == PTDB's E1 F_MAX=4000) so dataset must be part of
# the key -- this bit us once already in compile_best_results.py/
# analyze_512bins_inference.py as an N_Bins-only version of the same mistake.
#   E1: full repr up to fs/2   E2: full repr up to fs/4
#   E3: upper-frequency truncation at fs/4   E4: upper-bin truncation at fs/4
TECHNIQUE_MAP = {
    # PTDB, 8kHz native
    ("ptdb", "Original-linear-bins", 4000): "E1", ("ptdb", "Mel-aligned-bins", 4000): "E1",
    ("ptdb", "Original-linear-bins", 2000): "E2", ("ptdb", "Mel-aligned-bins", 2000): "E2",
    ("ptdb", "Original-linear-upper-frequency-truncated", 2000): "E3",
    ("ptdb", "Mel-aligned-upper-frequency-truncated", 2000): "E3",
    ("ptdb", "Original-linear-upper-bins-truncated", 2000): "E4",
    ("ptdb", "Mel-aligned-upper-bins-truncated", 2000): "E4",
    # MIR-1K, 16kHz native
    ("mir1k", "Original-linear-bins", 8000): "E1", ("mir1k", "Mel-aligned-bins", 8000): "E1",
    ("mir1k", "Original-linear-bins", 4000): "E2", ("mir1k", "Mel-aligned-bins", 4000): "E2",
    ("mir1k", "Original-linear-upper-frequency-truncated", 4000): "E3",
    ("mir1k", "Mel-aligned-upper-frequency-truncated", 4000): "E3",
    ("mir1k", "Original-linear-upper-bins-truncated", 4000): "E4",
    ("mir1k", "Mel-aligned-upper-bins-truncated", 4000): "E4",
}


def parse_sr_fl_bins(name):
    m = _SR_FL_BINS_NEW.fullmatch(name)
    if m:
        return int(m.group(1)), int(m.group(2)), int(m.group(3))
    m = _SR_FL_BINS_OLD.fullmatch(name)
    if m:
        return int(m.group(1)), int(m.group(2)), int(m.group(3))
    return None


def get_last_metrics(result_file):
    last = None
    with open(result_file) as f:
        for line in f:
            nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", line)
            if len(nums) < 6:
                continue
            try:
                last = {
                    "iteration": int(float(nums[0])),
                    "RPA": float(nums[1]),
                    "RCA": float(nums[2]),
                    "OA": float(nums[3]),
                    "VR": float(nums[4]),
                    "VFA": float(nums[5]),
                }
            except (ValueError, TypeError):
                continue
    return last


def get_test_metrics(run_dir):
    """Extract clean condition for ranking, plus all noise conditions for reference."""
    test_eval = run_dir / "test_eval.csv"
    if not test_eval.exists():
        return {}

    metrics = {}
    conditions_found = {}

    with open(test_eval, newline="") as f:
        for row in csv.DictReader(f):
            cond = row["condition"]
            conditions_found[cond] = {
                "RPA": round(float(row["RPA"]), 4),
                "RCA": round(float(row["RCA"]), 4),
                "OA": round(float(row["OA"]), 4),
                "VR": round(float(row["VR"]), 4),
                "VFA": round(float(row["VFA"]), 4),
            }

    # Primary ranking metric: clean condition
    if "clean" in conditions_found:
        c = conditions_found["clean"]
        metrics.update({
            "test_RPA": c["RPA"],
            "test_RCA": c["RCA"],
            "test_OA": c["OA"],
            "test_VR": c["VR"],
            "test_VFA": c["VFA"],
        })

    # Also capture noise condition RPA for visibility
    for snr in ["snr00", "snr05", "snr10", "snr20"]:
        if snr in conditions_found:
            metrics[f"test_RPA_{snr}"] = conditions_found[snr]["RPA"]

    return metrics


def main():
    rows = []
    skipped = 0

    for result_file in sorted(RUNS_DIR.rglob("result.txt")):
        run_dir = result_file.parent
        parts = run_dir.parts

        # Fixed layout relative to RUNS_DIR (runs_seeded/pitch_th_X/...), set
        # by train_single.py's logdir construction -- reading by position
        # instead of searching for known literal strings (the old "clean"/
        # "clean_noisy"/"ptdb_sr" checks) is what lets this generalize to
        # any dataset/datatype without hardcoding, e.g. MIR-1K's
        # "polyphonic"/"clean_vocals" datatypes and "mir1k_sr16000" tag.
        #   .../pitch_th_X/{dataset}_sr{rate}/{datatype}/{fmax_tag}/{feature_dir}/{pitch_dir}
        #
        # `datatype` itself is now a variable number of segments -- e.g.
        # plain "clean" (PTDB), or "clean+noisy/polyphonic" / "clean/polyphonic"
        # (MIR-1K, now that a train-condition level sits above the poly/
        # clean_vocals split -- see MIR1K_CLEAN_NOISY_CAMPAIGN_STATUS.md) --
        # so the dataset_tag segment is found by scanning rather than a
        # fixed parts[-5] offset.
        if len(parts) < 5:
            skipped += 1
            continue

        pitch_dir = parts[-1]
        feature_dir = parts[-2]
        sr_fl_folder = parts[-3]

        dataset_tag_idx = None
        for i in range(len(parts) - 4, -1, -1):
            if _DATASET_SR.match(parts[i]):
                dataset_tag_idx = i
                break
        if dataset_tag_idx is None:
            skipped += 1
            continue

        m = _DATASET_SR.match(parts[dataset_tag_idx])
        dataset, SR = m.group(1), int(m.group(2))
        # e.g. ("clean",) -> "clean" ; ("clean+noisy","polyphonic") -> "clean+noisy_polyphonic"
        data_condition = "_".join(parts[dataset_tag_idx + 1:-3])

        m = _FEATURE_DIR.fullmatch(feature_dir)
        if not m:
            skipped += 1
            continue
        feature_type, lr, batch_size = m.group(1), float(m.group(2)), int(m.group(3))
        seed_tag = m.group(4) or ""  # "" = default seed=42, untagged (see _FEATURE_DIR's comment)

        parsed = parse_sr_fl_bins(sr_fl_folder)
        if parsed is None:
            skipped += 1
            continue
        f_max, window_length, n_bins = parsed

        window_ms = round(window_length / SR * 1000)

        pm = _PITCH_DIR.fullmatch(pitch_dir)
        if not pm:
            skipped += 1
            continue
        alpha, gamma = float(pm.group(1)), float(pm.group(2))

        metrics = get_last_metrics(result_file)
        if metrics is None:
            skipped += 1
            continue

        test_metrics = get_test_metrics(run_dir)

        rows.append({
            # meta (sweep-constant, kept per-row here only to derive the
            # filename below, then dropped before writing the CSV)
            "_hop_ms": HOP_MS_BY_DATASET.get(dataset, 10),
            "_dataset": dataset,
            "_data_condition": data_condition,
            "_SR": SR,
            "_batch_size": batch_size,
            "_window_ms": window_ms,
            # per-row columns
            "FL": window_length,
            "N_Bins": n_bins,
            "F_MAX": f_max,
            "feature_type": feature_type,
            "technique": TECHNIQUE_MAP.get((dataset, feature_type, f_max), ""),
            "batch_size": batch_size,
            "seed_tag": seed_tag,
            "alpha": alpha,
            "gamma": gamma,
            **metrics,
            **test_metrics,
        })

    print(f"[{datetime.now().isoformat()}] Rows: {len(rows)}  Skipped: {skipped}")
    if not rows:
        print("No rows collected yet -- nothing to write.")
        return

    # Group by (dataset, data_condition, window_ms, batch_size, N_Bins) --
    # e.g. ptdb/clean, mir1k/polyphonic, mir1k/clean_vocals each get their
    # OWN csv, and now window_ms/batch_size/N_Bins do too (2026-09-19,
    # for the PTDB batch-size/bin-count matrix -- these files
    # feed plotting/publishing directly, so "one file = one exact
    # (wl, bs, n_bins, condition) combination" matters more than it used
    # to). Merging different combinations into one file would silently
    # average/rank across incompatible sweeps (different SR, hop, F_MAX
    # scale, even different datasets).
    groups = {}
    for r in rows:
        groups.setdefault(
            (r["_dataset"], r["_data_condition"], int(r["_window_ms"]), int(r["_batch_size"]), int(r["N_Bins"])),
            [],
        ).append(r)

    for (dataset, data_condition, window_ms, batch_size, n_bins), group_rows in sorted(groups.items()):
        group_rows.sort(key=lambda r: (r["F_MAX"], r["feature_type"], r["alpha"]))

        hop_ms = group_rows[0]["_hop_ms"]
        SR = group_rows[0]["_SR"]

        stem = (f"{dataset}_{hop_ms}ms_sr{SR}_{data_condition}"
                f"_wl{window_ms}ms_bs{batch_size}_{n_bins}bins_validation_results")
        csv_path = OUT_DIR / f"{stem}.csv"

        clean_rows = [{k: r.get(k, "") for k in COLUMNS} for r in group_rows]
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=COLUMNS)
            writer.writeheader()
            writer.writerows(clean_rows)

        print(f"Saved CSV: {csv_path}  ({len(group_rows)} rows)")


if __name__ == "__main__":
    main()

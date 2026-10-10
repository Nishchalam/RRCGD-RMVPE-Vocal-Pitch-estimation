"""
Final results compilation for the clean-data sweep(s), across every
dataset/data-condition export_validation_csv.py has produced.

Reads each *_validation_results CSV (one row per completed (experiment,
alpha, n_bins), with both validation and test-set metrics already in it)
instead of re-walking runs_seeded/ -- much cheaper once the CSV already
exists, and it's kept fresh automatically every 10 completed jobs by
auto_export_every_10.sh.

IMPORTANT: results are grouped by (feature_type, F_MAX, N_Bins, batch_size)
-- NOT just (feature_type, F_MAX) -- since a single CSV can contain rows
from multiple n_bins sweeps (e.g. the 512-bin production sweep AND a
256-bin batch-size ablation) mixed together. Grouping without N_Bins
silently merges e.g. 23 alphas of 512-bin "Original-linear-bins,F_MAX=4000"
with 15 alphas of 256-bin "Original-linear-bins,F_MAX=4000" into one
inflated 38-alpha row. batch_size is included in the group key too (even
though export_validation_csv.py already files different batch sizes into
separate CSVs via a data_condition tag) as defensive redundancy -- see the
2026-09-19 incident where the 256-bin sweep silently ran at bs=32 instead
of the finalized bs=8 for ~155 jobs before anyone noticed.

Dataset is ALSO part of the technique-mapping key, not just F_MAX: MIR-1K's
E2 (F_MAX=4000, its fs/4 at 16kHz) has the same numeric F_MAX as PTDB's E1
(F_MAX=4000, its fs/2 at 8kHz) -- a second instance of the same class of
bug. Since export_validation_csv.py now writes one CSV per (dataset,
data_condition) already, dataset is inferred from each CSV's filename
prefix rather than needing a redundant per-row column.

By default processes EVERY *_validation_results.csv found (not just one --
autodetecting only the alphabetically-last match would silently ignore
the others now that MIR-1K and PTDB coexist), writing one
best_results_<dataset>_<data_condition>_wl<wl>ms_bs<bs>_<n_bins>bins.csv
per (dataset, data_condition, window_ms, batch_size, n_bins) combination
found -- window_ms/batch_size/n_bins are read straight from each input
filename (export_validation_csv.py, 2026-09-19: always-explicit naming,
one exact combination per file) rather than re-derived from row data,
though N_Bins is still cross-checked against the data as a safety net.
"""
import argparse
import csv
import re
from pathlib import Path

METRICS = ["RPA", "RCA", "OA", "VR", "VFA"]

# E4-linear ("Original-linear-upper-bins-truncated") retired 2026-09-19:
# on a LINEAR grid, "keep bins with center frequency <= fs/4" (E3-linear)
# and "keep the first half of bins by index" (E4-linear) select the exact
# same bin set -- verified byte-identical on both datasets (PTDB:
# 139/139 paired (config,alpha) points; MIR-1K: 46/46 for both
# clean_polyphonic and clean_clean_vocals; max|diff|=0.000000 in every
# case). NOT true for mel-aligned grids (E4-mel is genuinely different
# from E3-mel and stays). See ptdb_10ms_wl32ms_plan.md for the analysis.
#
# Rows with this feature_type are dropped in compile_one_csv() before
# grouping, for every dataset -- for PTDB this is a no-op (the duplicate
# run directories were deleted, so no such rows exist anymore); for
# MIR-1K it's the actual mechanism, since that sweep's checkpoints are
# gone and its raw *_validation_results.csv (92 duplicate rows) can't be
# regenerated -- only the compiled best_results_*.csv drops the row, the
# raw CSV stays as committed audit trail.
RETIRED_DUPLICATE_FEATURE_TYPES = {"Original-linear-upper-bins-truncated"}

# (dataset, feature_type, F_MAX) -> (technique, bin_family) -- mirrors each
# dataset's own experiments_*.yaml representations.
TECHNIQUE_MAP = {
    # PTDB, 8kHz native
    ("ptdb", "Original-linear-bins", 4000): ("E1", "linear"),
    ("ptdb", "Mel-aligned-bins", 4000): ("E1", "mel"),
    ("ptdb", "Original-linear-bins", 2000): ("E2", "linear"),
    ("ptdb", "Mel-aligned-bins", 2000): ("E2", "mel"),
    ("ptdb", "Original-linear-upper-frequency-truncated", 2000): ("E3", "linear"),
    ("ptdb", "Mel-aligned-upper-frequency-truncated", 2000): ("E3", "mel"),
    ("ptdb", "Mel-aligned-upper-bins-truncated", 2000): ("E4", "mel"),
    # MIR-1K, 16kHz native
    ("mir1k", "Original-linear-bins", 8000): ("E1", "linear"),
    ("mir1k", "Mel-aligned-bins", 8000): ("E1", "mel"),
    ("mir1k", "Original-linear-bins", 4000): ("E2", "linear"),
    ("mir1k", "Mel-aligned-bins", 4000): ("E2", "mel"),
    ("mir1k", "Original-linear-upper-frequency-truncated", 4000): ("E3", "linear"),
    ("mir1k", "Mel-aligned-upper-frequency-truncated", 4000): ("E3", "mel"),
    ("mir1k", "Mel-aligned-upper-bins-truncated", 4000): ("E4", "mel"),
}

# export_validation_csv.py (2026-09-19: always-explicit naming) writes
# <dataset>_<hop>ms_sr<rate>_<condition>_wl<wl>ms_bs<bs>_<n>bins_validation_results.csv --
# every file is now homogeneous in window-length/batch-size/bin-count by
# construction (one exact combination per file), so those three values
# come straight from the filename instead of needing re-derivation from
# the row data.
_CSV_NAME = re.compile(
    r"^([a-zA-Z0-9]+)_\d+ms_sr\d+_([a-zA-Z0-9_+]+)_wl(\d+)ms_bs(\d+)_(\d+)bins_validation_results\.csv$"
)


def parse_validation_csv_filename(csv_path: Path):
    """Recover (dataset, data_condition, window_ms, batch_size, n_bins) from
    export_validation_csv.py's naming convention, e.g.
    'ptdb_10ms_sr8000_clean_wl32ms_bs8_512bins_validation_results.csv'
    -> ('ptdb', 'clean', 32, 8, 512)."""
    m = _CSV_NAME.match(csv_path.name)
    if not m:
        raise SystemExit(
            f"Can't parse {csv_path.name!r} -- expected export_validation_csv.py's "
            f"<dataset>_<hop>ms_sr<rate>_<condition>_wl<wl>ms_bs<bs>_<n>bins_validation_results.csv"
        )
    dataset, condition, wl, bs, n_bins = m.groups()
    return dataset, condition, int(wl), int(bs), int(n_bins)


def find_csvs(base: Path, explicit: str | None) -> list[Path]:
    if explicit:
        return [base / explicit]
    matches = sorted(base.glob("*_validation_results.csv"))
    if not matches:
        raise SystemExit(f"No *_validation_results.csv found in {base} -- run export_validation_csv.py first.")
    return matches


def compile_for_n_bins(rows_for_n_bins, dataset, data_condition, n_bins, window_ms, batch_size, base):
    """Compile and write best_results_<dataset>_<data_condition>_wl<wl>ms_bs<bs>_<n_bins>bins.csv."""
    groups = {}  # (feature_type, F_MAX, batch_size) -> list of row dicts
    for row in rows_for_n_bins:
        key = (row["feature_type"], int(row["F_MAX"]), int(row["batch_size"]))
        groups.setdefault(key, []).append(row)

    out_rows = []
    print(f"\n=== {dataset}/{data_condition}  N_Bins={n_bins} ===")
    for key, rows in sorted(groups.items()):
        technique, bin_family = TECHNIQUE_MAP.get((dataset, key[0], key[1]), ("?", "?"))
        n_total = len(rows)

        def rpa_for(r):
            return float(r["test_RPA"]) if r.get("test_RPA") else float(r["RPA"])

        def source_for(r):
            return "test" if r.get("test_RPA") else "val (test pending)"

        best = max(rows, key=rpa_for)
        n_test_done = sum(1 for r in rows if r.get("test_RPA"))

        label = f"{technique}_{bin_family}"
        print(f"{label:10s}  {feature_type_label(key)[:45]:45s}  "
              f"{n_test_done}/{n_total} test-evaluated  "
              f"best alpha={best['alpha']:<5}  RPA={rpa_for(best):.4f}  source={source_for(best)}")

        out_rows.append({
            "technique": technique,
            "bin_family": bin_family,
            "feature_type": key[0],
            "F_MAX": key[1],
            "N_Bins": n_bins,
            "batch_size": key[2],
            "n_alphas_evaluated": n_test_done,
            "n_alphas_total": n_total,
            "best_alpha": best["alpha"],
            "RPA": rpa_for(best),
            "RCA": float(best["test_RCA"]) if best.get("test_RCA") else float(best["RCA"]),
            "OA": float(best["test_OA"]) if best.get("test_OA") else float(best["OA"]),
            "VR": float(best["test_VR"]) if best.get("test_VR") else float(best["VR"]),
            "VFA": float(best["test_VFA"]) if best.get("test_VFA") else float(best["VFA"]),
            "source": source_for(best),
        })

    if not out_rows:
        print(f"No rows found for {dataset}/{data_condition} N_Bins={n_bins}.")
        return None

    out_path = base / f"best_results_{dataset}_{data_condition}_wl{window_ms}ms_bs{batch_size}_{n_bins}bins.csv"
    with open(out_path, "w", newline="") as f:
        fieldnames = ["technique", "bin_family", "feature_type", "F_MAX", "N_Bins", "batch_size",
                      "n_alphas_evaluated", "n_alphas_total", "best_alpha"] + METRICS + ["source"]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(out_rows)

    print(f"Wrote {len(out_rows)} experiment(s) to {out_path}")
    return out_path


def compile_one_csv(csv_path: Path, base: Path, only_n_bins: int | None):
    dataset, data_condition, file_wl, file_bs, file_n_bins = parse_validation_csv_filename(csv_path)
    print(f"\nReading {csv_path}  (dataset={dataset}, data_condition={data_condition}, "
          f"wl={file_wl}ms, bs={file_bs}, n_bins={file_n_bins})")

    if only_n_bins is not None and only_n_bins != file_n_bins:
        print(f"  Skipping -- this file is N_Bins={file_n_bins}, not the requested {only_n_bins}.")
        return

    # Each file is homogeneous in (window_ms, batch_size, N_Bins) by
    # construction (export_validation_csv.py's naming), but re-derive
    # N_Bins from the row data too as a defensive check -- if a stray
    # mixed-N_Bins file ever showed up (e.g. from an older export before
    # this convention), better to bucket correctly than silently trust
    # a filename that turned out to lie.
    by_n_bins = {}
    n_dropped = 0
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            if row["feature_type"] in RETIRED_DUPLICATE_FEATURE_TYPES:
                n_dropped += 1
                continue
            by_n_bins.setdefault(int(row["N_Bins"]), []).append(row)

    if n_dropped:
        print(f"  Dropped {n_dropped} row(s) of retired duplicate technique(s) "
              f"{sorted(RETIRED_DUPLICATE_FEATURE_TYPES)} (see RETIRED_DUPLICATE_FEATURE_TYPES).")

    if set(by_n_bins) != {file_n_bins}:
        print(f"  WARNING: filename says N_Bins={file_n_bins} but data contains {sorted(by_n_bins)} -- "
              f"compiling each found N_Bins separately rather than trusting the filename.")

    for n_bins, rows in sorted(by_n_bins.items()):
        compile_for_n_bins(rows, dataset, data_condition, n_bins, file_wl, file_bs, base)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=None,
                     help="Explicit validation_results CSV filename (default: process every "
                          "*_validation_results.csv found)")
    ap.add_argument("--n-bins", type=int, default=None,
                     help="Only compile for this N_Bins value (default: compile one file per N_Bins found)")
    args = ap.parse_args()

    base = Path(__file__).parent
    csv_paths = find_csvs(base, args.csv)

    for csv_path in csv_paths:
        compile_one_csv(csv_path, base, args.n_bins)


def feature_type_label(key):
    return f"{key[0]} (F_MAX={key[1]})"


if __name__ == "__main__":
    main()

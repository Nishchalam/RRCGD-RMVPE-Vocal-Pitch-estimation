#!/usr/bin/env python3
"""Per-SNR-condition breakdown for each technique's current best-alpha
checkpoint, at both 256 and 512 bins -- no re-evaluation, just reads the
test_eval.csv each checkpoint already produced (condition rows: clean,
snr00/05/10/20, overall). Reads whichever best_results_ptdb_clean_noisy_
wl32ms_bs8_{256,512}bins.csv compile_best_results.py last wrote, so run
that first (or let the scheduler, which already does).

Output: ptdb_clean_noisy_wl32ms_bs8_results_by_condition.csv -- one row
per (technique, N_Bins, condition). The best_alpha per technique will
change as the clean+noisy sweep fills in, so re-run this after any
compile_best_results.py refresh to keep it current.
"""
import csv
from pathlib import Path

BASE = Path(__file__).resolve().parent
RUNS = BASE / "runs_seeded/pitch_th_0.5/ptdb_sr8000/clean_noisy"
BEST_CSVS = ["best_results_ptdb_clean_noisy_wl32ms_bs8_256bins.csv",
             "best_results_ptdb_clean_noisy_wl32ms_bs8_512bins.csv"]
OUT = BASE / "ptdb_clean_noisy_wl32ms_bs8_results_by_condition.csv"


def run_dir_for(feature_type, f_max, n_bins, batch_size, alpha):
    return RUNS / f"{f_max}_FL_256_bins{n_bins}" / f"{feature_type}_0.001_{batch_size}" / f"Pitch_FL_{alpha}_0.0"


def main():
    out_rows = []
    for name in BEST_CSVS:
        p = BASE / name
        if not p.exists():
            print(f"skip {name}: not found (run compile_best_results.py first)")
            continue
        for r in csv.DictReader(open(p)):
            rd = run_dir_for(r["feature_type"], r["F_MAX"], r["N_Bins"], r["batch_size"], r["best_alpha"])
            te = rd / "test_eval.csv"
            if not te.exists():
                print(f"  WARNING: no test_eval.csv at {rd} (best_results.csv is stale?) -- skipping")
                continue
            for row in csv.DictReader(open(te)):
                out_rows.append({
                    "technique": r["technique"], "bin_family": r["bin_family"],
                    "feature_type": r["feature_type"], "F_MAX": r["F_MAX"], "N_Bins": r["N_Bins"],
                    "best_alpha": r["best_alpha"],
                    **row,
                })

    if not out_rows:
        print("No rows -- nothing written.")
        return

    fields = ["technique", "bin_family", "feature_type", "F_MAX", "N_Bins", "best_alpha",
              "condition", "n_files", "n_frames", "RPA", "RCA", "OA", "VR", "VFA"]
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(out_rows)
    print(f"Wrote {len(out_rows)} rows to {OUT}")


if __name__ == "__main__":
    main()

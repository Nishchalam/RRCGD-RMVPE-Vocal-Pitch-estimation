#!/usr/bin/env python3
"""
Evaluate one saved optimized RRCGD-RMVPE checkpoint on held-out PTDB test data.

Reports:
  clean
  every SNR condition found in the test split
  overall = all test conditions combined, frame-weighted

The model checkpoint is never selected using this test evaluation.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
from mir_eval.melody import (
    raw_pitch_accuracy,
    raw_chroma_accuracy,
    overall_accuracy,
    voicing_recall,
    voicing_false_alarm,
    to_cent_voicing,
)

from src import PTDB, to_local_average_cents


def evaluate_item(model, item, hop_ms, device):
    roots = item["roots"].to(device)
    labels = item["pitch"].to(device)

    if roots.dim() == 2:
        roots = roots.unsqueeze(0)

    with torch.no_grad():
        logits = model(roots).squeeze(0)
        pred = torch.sigmoid(logits)

    T = min(labels.shape[0], pred.shape[0])
    pred = pred[:T].cpu().numpy()
    labels = labels[:T].cpu().numpy()

    cents_p = to_local_average_cents(pred, None, 0.5)
    cents_y = to_local_average_cents(labels, None, 0.5)

    est = np.array(
        [10 * 2 ** (c / 1200) if c else 0 for c in cents_p],
        dtype=float,
    )
    ref = np.array(
        [10 * 2 ** (c / 1200) if c else 0 for c in cents_y],
        dtype=float,
    )

    times = np.arange(T, dtype=float) * hop_ms / 1000.0

    return to_cent_voicing(times, ref, times, est)


def metrics(rv, rc, ev, ec):
    return {
        "RPA": float(raw_pitch_accuracy(rv, rc, ev, ec)),
        "RCA": float(raw_chroma_accuracy(rv, rc, ev, ec)),
        "OA": float(overall_accuracy(rv, rc, ev, ec)),
        "VR": float(voicing_recall(rv, ev)),
        "VFA": float(voicing_false_alarm(rv, ev)),
    }


def evaluate(model_path, data_path, hop_ms, device):
    model = torch.load(
        model_path,
        map_location=device,
        weights_only=False,
    ).to(device).eval()

    ds = PTDB(
        str(data_path),
        hop_length=hop_ms,
        sequence_length=None,
        groups=["test"],
        use_cache=True,
    )

    grouped = defaultdict(lambda: [[], [], [], []])
    all_arrays = [[], [], [], []]
    file_counts = defaultdict(int)

    for i in range(len(ds)):
        rv, rc, ev, ec = evaluate_item(model, ds[i], hop_ms, device)

        condition = ds.condition_at(i)
        grouped[condition][0].append(rv)
        grouped[condition][1].append(rc)
        grouped[condition][2].append(ev)
        grouped[condition][3].append(ec)
        file_counts[condition] += 1

        for j, arr in enumerate((rv, rc, ev, ec)):
            all_arrays[j].append(arr)

    rows = []

    def make_row(condition, arrays, n_files):
        arrays = [np.concatenate(x) for x in arrays]
        row = {
            "condition": condition,
            "n_files": n_files,
            "n_frames": len(arrays[0]),
        }
        row.update(metrics(*arrays))
        return row

    for condition in sorted(grouped):
        rows.append(
            make_row(
                condition,
                grouped[condition],
                file_counts[condition],
            )
        )

    # Overall is computed from all test frames, not as the mean
    # of the clean/SNR percentages.
    rows.append(
        make_row(
            "overall",
            all_arrays,
            len(ds),
        )
    )

    del model
    if torch.cuda.is_available():
        torch.cuda.empty_cache()

    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--data", required=True)
    ap.add_argument("--hop", required=True, type=int)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    device = torch.device(
        "cuda" if torch.cuda.is_available() else "cpu"
    )

    print("=" * 80)
    print(f"MODEL : {args.model}")
    print(f"DATA  : {args.data}")
    print(f"HOP   : {args.hop} ms")
    print(f"DEVICE: {device}")
    print("=" * 80)

    rows = evaluate(
        Path(args.model),
        Path(args.data),
        args.hop,
        device,
    )

    fields = [
        "condition",
        "n_files",
        "n_frames",
        "RPA",
        "RCA",
        "OA",
        "VR",
        "VFA",
    ]

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"{row['condition']:>8s} "
            f"n={row['n_files']:4d} "
            f"frames={row['n_frames']:8d} "
            + " ".join(
                f"{k}={row[k]:.4f}"
                for k in ("RPA", "RCA", "OA", "VR", "VFA")
            )
        )

    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()

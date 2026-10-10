#!/usr/bin/env python3
"""
Environment/determinism fingerprint for cross-machine reproducibility
debugging. Run this on BOTH machines (same gd_exps env activated) and
diff the two JSON reports -- any difference in the "environment" block is
a candidate explanation for non-matching results at the same
(seed, config, alpha); the "determinism_selftest" block checks whether
determinism holds AT ALL on this one machine, independent of the other.

Usage:
    conda activate gd_exps   # or however this env is activated on that box
    python3 check_repro_env.py
    # writes env_report_<hostname>.json in the current directory

Then compare env_report_<machine_a>.json vs env_report_<machine_b>.json --
diff the whole file (a plain `diff` or `jq` compare is enough) rather than
eyeballing selected fields, since it's easy to overlook a silently
different transitive dependency (e.g. a numpy or cuDNN version pulled in
differently by conda's solver on different platforms/GPU generations).
"""
import json
import platform
import socket
import subprocess
import sys
from pathlib import Path


def run(cmd):
    try:
        return subprocess.run(
            cmd, shell=True, capture_output=True, text=True, timeout=15
        ).stdout.strip()
    except Exception as e:
        return f"[error: {e}]"


def get_environment_info():
    info = {
        "hostname": socket.gethostname(),
        "python_version": sys.version,
        "platform": platform.platform(),
    }

    try:
        import torch
        info["torch_version"] = torch.__version__
        info["torch_cuda_build"] = torch.version.cuda
        info["cudnn_version"] = torch.backends.cudnn.version()
        info["cuda_available"] = torch.cuda.is_available()
        if torch.cuda.is_available():
            info["gpu_name"] = torch.cuda.get_device_name(0)
            info["gpu_compute_capability"] = list(torch.cuda.get_device_capability(0))
            info["gpu_count"] = torch.cuda.device_count()
        info["cudnn_deterministic_flag"] = torch.backends.cudnn.deterministic
        info["cudnn_benchmark_flag"] = torch.backends.cudnn.benchmark
        info["deterministic_algorithms_enabled"] = torch.are_deterministic_algorithms_enabled()
    except ImportError as e:
        info["torch_import_error"] = str(e)

    try:
        import numpy
        info["numpy_version"] = numpy.__version__
    except ImportError:
        info["numpy_version"] = None

    try:
        import scipy
        info["scipy_version"] = scipy.__version__
    except ImportError:
        info["scipy_version"] = None

    info["nvidia_driver_version"] = run(
        "nvidia-smi --query-gpu=driver_version --format=csv,noheader -i 0"
    )
    info["nvcc_version"] = run("nvcc --version 2>/dev/null | grep release") or "[nvcc not found -- using conda's cudatoolkit, not a system install]"
    info["conda_env_name"] = run("echo $CONDA_DEFAULT_ENV")
    info["full_pip_freeze"] = run("pip freeze 2>/dev/null")

    return info


def determinism_selftest():
    """Confirms determinism holds on THIS machine alone -- run the exact
    same op twice with the same seed and check bit-exact agreement. This
    is a prerequisite for cross-machine matching but doesn't by itself
    prove cross-machine matching (see module docstring)."""
    import torch

    def run_once():
        torch.manual_seed(42)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(42)
        torch.use_deterministic_algorithms(True, warn_only=True)
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        # A tiny conv + matmul stack -- representative of the actual model's
        # building blocks (E2E uses Conv2d + GRU + Linear), not just a raw
        # matmul, since conv algorithm selection is a common source of
        # cross-GPU-generation nondeterminism.
        x = torch.randn(4, 3, 32, 32, device=device)
        conv = torch.nn.Conv2d(3, 8, 3, padding=1).to(device)
        y = conv(x)
        z = torch.matmul(y.flatten(1), y.flatten(1).T)
        return z.detach().cpu()

    a = run_once()
    b = run_once()
    return {
        "same_machine_bit_exact": torch.equal(a, b),
        "max_abs_diff_if_any": float((a - b).abs().max()) if not torch.equal(a, b) else 0.0,
        "sample_values_run1": a[0, :3].tolist(),
        "sample_values_run2": b[0, :3].tolist(),
    }


def main():
    report = {
        "environment": get_environment_info(),
        "determinism_selftest": determinism_selftest(),
    }

    hostname = report["environment"]["hostname"]
    out_path = Path(f"env_report_{hostname}.json")
    out_path.write_text(json.dumps(report, indent=2))

    print(f"=== Environment report for {hostname} ===")
    for k, v in report["environment"].items():
        if k == "full_pip_freeze":
            continue  # too long for the console summary, still in the JSON
        print(f"  {k}: {v}")
    print()
    print("=== Determinism self-test (same machine, run twice) ===")
    for k, v in report["determinism_selftest"].items():
        print(f"  {k}: {v}")
    print()
    print(f"Full report (including pip freeze) written to: {out_path}")
    print(f"Copy this file from both machines and diff them -- e.g.:")
    print(f"  diff <(python3 -m json.tool env_report_MACHINE_A.json) <(python3 -m json.tool env_report_MACHINE_B.json)")


if __name__ == "__main__":
    main()

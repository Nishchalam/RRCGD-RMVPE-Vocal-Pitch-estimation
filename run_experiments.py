"""
run_experiments.py  (ACS of RRCGD)
────────────────────────────────────
Sequential experiment runner.  Reads experiments.yaml, expands each entry
over its alpha_list, and runs every (experiment × alpha) combination one
after another by spawning a fresh subprocess for each job.

Each subprocess inherits the current environment plus the EXP_* variables
that configure src/constants.py (SAMPLE_RATE, F_MAX, WINDOW_LENGTH,
BIN_METHOD) and train_single.py.

Usage
─────
    # Run all experiments
    python run_experiments.py

    # Dry-run: print what would be launched without executing anything
    python run_experiments.py --dry-run

    # Run specific experiments (comma-separated names)
    python run_experiments.py --only ptdb_fmax4000_wl2048_linear,poly_fmax8000_wl2048_linear

    # Use an alternative config file
    python run_experiments.py --config my_grid.yaml

    # Skip any job whose log directory already exists (resume after crash)
    python run_experiments.py --skip-existing
"""

import argparse
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import yaml   # pip install pyyaml

EVAL_SCRIPTS = {
    "ptdb":  Path(__file__).parent / "evaluate_ptdb_checkpoint.py",
    "mir1k": Path(__file__).parent / "evaluate_mir1k_checkpoint.py",
}


# ── Helpers ───────────────────────────────────────────────────────────────────

def hms(seconds: float) -> str:
    h = int(seconds // 3600)
    m = int((seconds % 3600) // 60)
    s = int(seconds % 60)
    return f"{h:02d}:{m:02d}:{s:02d}"

def resolve_f_max(f_max, sample_rate):
    """
    Resolve YAML f_max specification into an integer frequency in Hz.

    Supported forms:
        fs/2
        fs/4
        numeric values such as 4000
        numeric strings such as "4000"
    """
    if isinstance(f_max, (int, float)):
        return int(f_max)

    if isinstance(f_max, str):
        value = f_max.strip().lower()

        if value == "fs/2":
            return int(sample_rate // 2)

        if value == "fs/4":
            return int(sample_rate // 4)

        try:
            return int(float(value))
        except ValueError:
            raise ValueError(
                f"Unsupported f_max value: {f_max!r}. "
                f"Use 'fs/2', 'fs/4', or a numeric frequency in Hz."
            )

    raise ValueError(
        f"Unsupported f_max type: {type(f_max).__name__}"
    )

_METHOD_TAGS = {
    "mel_aligned": "Mel-aligned-bins",
    "linear_fmax": "Original-linear-bins",

    "mel_aligned_upper_frequency_truncated":
        "Mel-aligned-upper-frequency-truncated",

    "linear_fmax_upper_frequency_truncated":
        "Original-linear-upper-frequency-truncated",

    "mel_aligned_upper_bins_truncated":
        "Mel-aligned-upper-bins-truncated",

    "linear_fmax_upper_bins_truncated":
        "Original-linear-upper-bins-truncated",

    "maxpool": "Maxpooled-linear-bins",
    "mel_maxpool": "Mel-Maxpooled-linear-bins",
    "mel_adaptive": "Mel-Adpative-linear-bins",
}


def logdir_for(exp: dict, alpha: float, defaults: dict) -> str:
    """Reproduce the same path logic as train_single.py."""
    # Cast to float like train_single.py does (EXP_LR/EXP_GAMMA go through
    # float(os.environ[...])) — otherwise YAML values like "1e-3" that PyYAML
    # leaves as strings (no decimal point) produce a different-looking path
    # ("_1e-3_32") than what train_single.py actually creates ("_0.001_32").
    lr       = float(exp.get("learning_rate", defaults.get("learning_rate", 1e-3)))
    bs       = exp.get("batch_size",    defaults.get("batch_size", 8))
    gamma    = float(exp.get("gamma",   defaults.get("gamma", 0)))
    pitch_th = exp.get("pitch_th",      defaults.get("pitch_th", 0.5))
    sr = exp.get("sample_rate",defaults.get("sample_rate", 16000))
    f_max = resolve_f_max(exp["f_max"], sr)
    seed = exp.get("seed", defaults.get("seed", 42))
    no_seed = exp.get("no_seed", defaults.get("no_seed", False))
    deterministic = exp.get("deterministic", defaults.get("deterministic", True))
    repeat = exp.get("repeat", defaults.get("repeat", 1))

    method_tag  = _METHOD_TAGS.get(exp["bin_method"], exp["bin_method"])
    # Include SR in tag so ptdb (8k) and mir1k (16k) log dirs never collide
    dataset_tag = f"{exp['dataset']}_sr{sr}"
    bins_tag = f"bins{exp.get('n_bins', defaults.get('n_bins', 256))}"
    # fold into fmax_tag or method_tag, e.g.:
    fmax_tag = f"{f_max}_FL_{exp['window_length']}_{bins_tag}"
    # fmax_tag    = f"{exp['f_max']}_FL_{exp['window_length']}"
    # Only tag non-default seeds -- see train_single.py's matching comment;
    # this function must stay byte-for-byte in sync with train_single.py's
    # own logdir construction or --skip-existing checks the wrong path.
    seed_tag = "" if seed == 42 else f"_seed{seed}"
    # no_seed wins over seed_tag (an unseeded run's SEED value is
    # meaningless -- see train_single.py's matching comment).
    if no_seed:
        seed_tag = "_noseed"
    # repeat only tags when non-default (2, 3, ...), same "only tag what's
    # not the common case" convention as seed_tag -- a repeat=1 run's path
    # is identical to every pre-existing run, so nothing already on disk
    # ever gets reinterpreted by adding this axis.
    repeat_tag = "" if repeat == 1 else f"_rep{repeat}"
    # deterministic=False routes to its own runs root rather than tagging
    # the leaf dir -- these are a different, non-bit-exact regime and
    # shouldn't be findable by the same --skip-existing / compile_best_results
    # scans that assume runs_seeded/ is reproducible. deterministic=True
    # (the default, and every run so far) is completely unaffected.
    runs_root = os.environ.get("EXP_RUNS_ROOT") or ("runs_seeded" if deterministic else "runs_deterministic_false")
    return (
        f"{runs_root}/pitch_th_{pitch_th}"
        # f"runs_claude_data split/pitch_th_{pitch_th}"
        f"/{dataset_tag}/{exp['datatype']}/{fmax_tag}"
        f"/{method_tag}_{lr}_{bs}{seed_tag}{repeat_tag}"
        f"/Pitch_FL_{alpha}_{gamma}"
    )


def build_env(exp: dict, alpha: float, defaults: dict,
              data_paths: dict) -> dict:
    """Build the full environment for one (experiment × alpha) subprocess."""
    env = os.environ.copy()
    # RRCGDSpectrogram._chunked_rrcgd_sum allocates several same-sized-but-
    # not-identical tensors per chunk, thousands of times over a run; the
    # default CUDA caching allocator fragments under that pattern until it
    # can't find a contiguous block big enough even with plenty of total
    # free memory ("reserved but unallocated" in an OOM message is the
    # tell). expandable_segments avoids that by growing one segment instead
    # of hunting for a new contiguous block each time.
    env["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"
    # Without this, train_single.py's plain print() calls (validation RPA,
    # "New best model saved", "Finished:") sit in Python's stdio buffer and
    # never reach the log until the buffer fills or the process exits --
    # only tqdm's progress bar appears live, since tqdm flushes explicitly.
    # That silently breaks any live tail/monitor of the log.
    env["PYTHONUNBUFFERED"] = "1"

    # ── Sample rate (must be set BEFORE constants.py is imported in subprocess)
    sr = exp.get(
        "sample_rate",
        defaults.get("sample_rate", 16000)
    )

    env["EXP_SAMPLE_RATE"] = str(sr)

    f_max = resolve_f_max(
        exp["f_max"],
        sr
    )

    env["EXP_F_MAX"]         = str(f_max)
    env["EXP_WINDOW_LENGTH"] = str(exp["window_length"])
    env["EXP_BIN_METHOD"]    = exp["bin_method"]
    env["EXP_N_DENSE"]       = str(exp.get("n_dense", 2048))
    env["EXP_N_BINS"]        = str(exp.get("n_bins", defaults.get("n_bins", 256)))
    pitch_th = exp.get("pitch_th", defaults.get("pitch_th", 0.5))
    env["EXP_PITCH_TH"]      = str(pitch_th)

    # ── Dataset ───────────────────────────────────────────────────────────────
    env["EXP_DATASET"]       = exp["dataset"]
    env["EXP_DATATYPE"]      = exp["datatype"]
    env["EXP_DATA_PATH"]     = data_paths[exp["dataset"]]
    env["EXP_TRAIN_SPLITS"]  = repr(exp["train_splits"])
    env["EXP_VALID_SPLITS"]  = repr(exp["valid_splits"])

    # ── Loss / optimiser ──────────────────────────────────────────────────────
    env["EXP_ALPHA"]          = str(alpha)
    env["EXP_GAMMA"]          = str(exp.get("gamma",         defaults.get("gamma",         0)))
    env["EXP_LR"]             = str(exp.get("learning_rate", defaults.get("learning_rate", 1e-3)))
    env["EXP_BATCH_SIZE"]     = str(exp.get("batch_size",    defaults.get("batch_size",    8)))
    env["EXP_SEQ_L"]          = str(exp.get("seq_l",         defaults.get("seq_l",         2.55)))
    env["EXP_HOP_LENGTH"]     = str(exp.get("hop_length",    defaults.get("hop_length",    20)))
    env["EXP_HOP_MS"]          = str(exp.get("hop_length",    defaults.get("hop_length",    20)))
    env["EXP_CLEAN_ONLY"]      = "1" if exp.get("clean_only", defaults.get("clean_only", False)) else "0"
    env["EXP_USE_CACHE"]       = "1" if exp.get("use_cache", defaults.get("use_cache", True)) else "0"
    env["EXP_CLIP_GRAD_NORM"] = str(exp.get("clip_grad_norm",defaults.get("clip_grad_norm",3)))
    env["EXP_PATIENCE_EPOCHS"] = str(exp.get("patience_epochs", defaults.get("patience_epochs", 4)))

    deterministic = exp.get("deterministic", defaults.get("deterministic", True))
    env["EXP_SEED"]          = str(exp.get("seed", defaults.get("seed", 42)))
    env["EXP_DETERMINISTIC"] = "1" if deterministic else "0"
    env["EXP_NO_SEED"]       = "1" if exp.get("no_seed", defaults.get("no_seed", False)) else "0"
    env["EXP_REPEAT"]        = str(exp.get("repeat", defaults.get("repeat", 1)))
    if deterministic:
        # cuBLAS reads this only when it creates its handle, so it has to be
        # in the child's environment from the start — setting it inside the
        # child after CUDA init would be too late to take effect.
        env["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"

    return env


# ── Main ──────────────────────────────────────────────────────────────────────

def eval_test_set(exp: dict, alpha: float, ldir: str, env: dict, defaults: dict) -> str:
    """Run evaluate_ptdb_checkpoint.py on the job's best checkpoint, writing
    test_eval.csv into its own logdir. Reuses the exact EXP_* env the training
    subprocess ran with, so the features fed at test time always match what
    the model was actually trained on (see run_test_eval_10ms.py's docstring
    for why that reuse matters -- this sidesteps that class of bug entirely
    by never reconstructing the env from the directory name).

    Returns a short status string for the progress log; never raises -- a
    failed eval shouldn't be mistaken for a failed training job.
    """
    dataset = exp.get("dataset", "ptdb")
    if dataset not in EVAL_SCRIPTS:
        return f"eval=skipped(no eval script for dataset={dataset})"
    eval_script = EVAL_SCRIPTS[dataset]

    run_dir = Path(ldir)
    ckpts = sorted(
        run_dir.glob("model-1-*.pt"),
        key=lambda p: int(p.stem.split("-")[-1]),
    )
    if not ckpts:
        return "eval=skipped(no checkpoint)"
    model = ckpts[-1]

    output = run_dir / "test_eval.csv"
    hop_ms = env.get("EXP_HOP_MS", str(defaults.get("hop_length", 10)))

    # Reusing the training env keeps model-matching features (bin method,
    # f_max, window_length, n_bins) correct, but EXP_CLEAN_ONLY is a
    # train/valid-only concept -- dataset.py's PTDB/MIR1K.files() applies it
    # to EVERY group including "test", so a clean-trained model's test-eval
    # would otherwise silently drop most of the test split down to clean
    # only. Test-set evaluation must always cover the full test split
    # regardless of what the model was trained on.
    eval_env = dict(env)
    eval_env["EXP_CLEAN_ONLY"] = "0"

    eval_args = [
        sys.executable, str(eval_script),
        "--model", str(model),
        "--data", env["EXP_DATA_PATH"],
        "--hop", hop_ms,
        "--output", str(output),
    ]

    if dataset == "mir1k":
        # MIR1K's test group name depends on which train group produced the
        # checkpoint (train_audio -> test_audio, train_clean_vocals ->
        # test_clean_vocals) -- PTDB only ever has one test group so it
        # doesn't need this.
        train_splits = exp.get("train_splits", [])
        train_group = train_splits[0] if train_splits else "train_audio"
        test_group = train_group.replace("train_", "test_", 1)
        eval_args += ["--group", test_group]

    result = subprocess.run(
        eval_args,
        env=eval_env,
        check=False,
    )
    if result.returncode != 0:
        return f"eval=FAILED(rc={result.returncode})"

    # Surface the headline number (overall RPA) in the progress log without
    # re-parsing the csv format elsewhere.
    try:
        import csv as _csv
        with open(output) as f:
            rows = {r["condition"]: r for r in _csv.DictReader(f)}
        rpa = float(rows["overall"]["RPA"])
        return f"eval=OK(test_RPA={rpa:.4f})"
    except Exception:
        return "eval=OK"


def main():
    parser = argparse.ArgumentParser(description="RRCGD-RMVPE sequential experiment runner")
    parser.add_argument(
        "--config", default="experiments.yaml",
        help="Path to YAML experiment config (default: experiments-mel.yaml)"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print jobs that would run without executing them"
    )
    parser.add_argument(
        "--only", default="",
        help="Comma-separated experiment names to run (others skipped)"
    )
    parser.add_argument(
        "--skip-existing", action="store_true",
        help="Skip any job whose log directory already exists"
    )
    parser.add_argument(
        "--use-cache", choices=["config", "true", "false"], default="config",
        help="Override every experiment's use_cache: setting (default: "
             "use whatever the config says, else on). 'false' bypasses "
             "the on-disk feature cache entirely -- raw audio/.mat "
             "parsed fresh every __getitem__ call."
    )
    parser.add_argument(
        "--no-eval-test", action="store_true",
        help="Skip the automatic post-training test-set evaluation "
             "(overrides eval_test:true in the config)"
    )
    parser.add_argument(
        "--deterministic", choices=["config", "true", "false"], default="config",
        help="Override every experiment's deterministic: setting (default: "
             "use whatever the config says). 'false' trades bit-exact repro "
             "for cuDNN's memory-efficient autotuned kernels -- useful when "
             "running multiple jobs concurrently on one GPU. seed: still "
             "applies either way, so init weights and data order stay fixed."
    )
    parser.add_argument(
        "--batch-size", type=int, default=None,
        help="Override every experiment's batch_size (default: use whatever "
             "the config says)."
    )
    parser.add_argument(
        "--patience-epochs", type=int, default=None,
        help="Override every experiment's early-stopping patience, in "
             "epochs (default: use whatever the config says, else 4)."
    )
    parser.add_argument(
        "--seed", type=int, default=None,
        help="Override every experiment's seed (default: use whatever the "
             "config says, else 42). Non-default seeds get a _seedN tag in "
             "the run dir path so they never collide with the default-seed "
             "run's checkpoint."
    )
    parser.add_argument(
        "--repeat", type=int, default=None,
        help="Override every experiment's repeat: index (default: use "
             "whatever the config says, else 1). Only tags the run dir "
             "when != 1 (_repN suffix) -- meant for re-running the exact "
             "same (seed, config) combo multiple times under "
             "deterministic=false to observe run-to-run variance from "
             "cuDNN/kernel nondeterminism alone."
    )
    parser.add_argument(
        "--no-seed", action="store_true",
        help="Skip set_seed() and the DataLoader's manual generator/"
             "worker_init_fn entirely -- plain unseeded PyTorch/NumPy/"
             "random RNG state, matching the PreRepro codebase "
             "(rrcgd_rmvpe_final.zip), which never seeded anything. "
             "Overrides --seed/--deterministic for the run itself (the "
             "env vars are still set, just unused by train_single.py)."
    )
    args = parser.parse_args()

    # ── Load config ───────────────────────────────────────────────────────────
    config_path = Path(args.config)
    if not config_path.exists():
        print(f"[ERROR] Config file not found: {config_path}")
        sys.exit(1)

    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    defaults    = cfg.get("defaults",    {})
    data_paths  = cfg.get("data_paths",  {})
    experiments = cfg.get("experiments", [])

    if args.deterministic != "config":
        defaults = dict(defaults)
        defaults["deterministic"] = (args.deterministic == "true")
        experiments = [
            {k: v for k, v in exp.items() if k != "deterministic"}
            for exp in experiments
        ]

    if args.batch_size is not None:
        defaults = dict(defaults)
        defaults["batch_size"] = args.batch_size
        experiments = [
            {k: v for k, v in exp.items() if k != "batch_size"}
            for exp in experiments
        ]

    if args.patience_epochs is not None:
        defaults = dict(defaults)
        defaults["patience_epochs"] = args.patience_epochs
        experiments = [
            {k: v for k, v in exp.items() if k != "patience_epochs"}
            for exp in experiments
        ]

    if args.seed is not None:
        defaults = dict(defaults)
        defaults["seed"] = args.seed
        experiments = [
            {k: v for k, v in exp.items() if k != "seed"}
            for exp in experiments
        ]

    if args.use_cache != "config":
        defaults = dict(defaults)
        defaults["use_cache"] = (args.use_cache == "true")
        experiments = [
            {k: v for k, v in exp.items() if k != "use_cache"}
            for exp in experiments
        ]

    if args.no_seed:
        defaults = dict(defaults)
        defaults["no_seed"] = True
        experiments = [
            {k: v for k, v in exp.items() if k != "no_seed"}
            for exp in experiments
        ]

    if args.repeat is not None:
        defaults = dict(defaults)
        defaults["repeat"] = args.repeat
        experiments = [
            {k: v for k, v in exp.items() if k != "repeat"}
            for exp in experiments
        ]

    only_set = {n.strip() for n in args.only.split(",") if n.strip()}

    # ── Build job list ────────────────────────────────────────────────────────
    jobs = []
    for exp in experiments:
        if only_set and exp["name"] not in only_set:
            continue
        # Per-experiment alpha_list takes priority over the global default
        alpha_list = exp.get("alpha_list", defaults.get("alpha_list", [6]))
        for alpha in alpha_list:
            jobs.append((exp, float(alpha)))

    total = len(jobs)
    print(f"\n{'━'*65}")
    print(f"  Config   : {config_path}")
    print(f"  Jobs     : {total}  ({len(experiments)} experiments × alphas)")
    print(f"  Mode     : {'DRY RUN' if args.dry_run else 'LIVE'}")
    print(f"  Started  : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'━'*65}\n")

    if total == 0:
        print("No jobs to run.  Check --only filter or experiments list.")
        return

    # ── Progress log ──────────────────────────────────────────────────────────
    # log_path = Path("experiment_run_paper_mel_aligned.log")
    log_path = Path("experiment_run_paper_linear.log")
    if not args.dry_run:
        with open(log_path, "a") as lf:
            lf.write(
                f"\n{'='*65}\n"
                f"Run started : {datetime.now().isoformat()}\n"
                f"Config      : {config_path}\n"
                f"Total jobs  : {total}\n"
            )

    # ── Execute ───────────────────────────────────────────────────────────────
    global_start = time.time()
    failed = []

    for job_idx, (exp, alpha) in enumerate(jobs, start=1):
        ldir  = logdir_for(exp, alpha, defaults)
        sr    = exp.get("sample_rate", defaults.get("sample_rate", 16000))
        label = (
            f"[{job_idx}/{total}]  {exp['name']}"
            f"  α={alpha}"
            f"  sr={sr}  fmax={exp['f_max']}  wl={exp['window_length']}"
            f"  bins={exp['bin_method']}"
            f"  hop={exp.get('hop_length', defaults.get('hop_length', 20))}ms"
            f"  data={'clean' if exp.get('clean_only', defaults.get('clean_only', False)) else 'clean+noisy'}"
        )

        # train_single.py's os.makedirs(logdir, exist_ok=True) runs before
        # any training happens, so the directory existing proves nothing --
        # a job that OOM'd or crashed on iteration 1 leaves an empty dir
        # behind too. A checkpoint only appears once a validation epoch
        # actually beat the running best, so its presence is what "this job
        # already succeeded" really means.
        if args.skip_existing and any(Path(ldir).glob("model-1-*.pt")):
            # Checkpoint exists, but has it actually been *evaluated*? A
            # training job that completed and then got interrupted (reboot,
            # manual stop, crash) before its post-training eval_test_set()
            # call finishes leaves exactly this state: real checkpoint, no
            # test_eval.csv. Treating that as "done" (the old behavior) is
            # a real, repeatedly-hit bug: a disk-based scheduler keeps
            # re-queuing the technique forever (it never reaches 23/23 on
            # disk), each relaunch SKIPs every alpha instantly and exits,
            # and the scheduler burns its launch slot on that GPU every
            # cycle without ever doing real work. Running just the missing
            # eval (not a full retrain) is cheap and closes the gap for
            # good instead of requiring a one-off manual fix each time.
            if not (Path(ldir) / "test_eval.csv").exists() and not args.no_eval_test \
                    and exp.get("eval_test", defaults.get("eval_test", True)):
                print(f"  SKIP   {label}  (checkpoint exists, no test_eval.csv -- evaluating only)")
                env = build_env(exp, alpha, defaults, data_paths)
                eval_status = eval_test_set(exp, alpha, ldir, env, defaults)
                print(f"  {eval_status}")
            else:
                print(f"  SKIP   {label}  (checkpoint exists)")
            continue

        print(f"\n{'─'*65}")
        print(f"  START  {label}")
        print(f"  logdir : {ldir}")
        print(f"{'─'*65}")

        if args.dry_run:
            env = build_env(exp, alpha, defaults, data_paths)
            print("  env overrides:")
            for k in sorted(k for k in env if k.startswith("EXP_")):
                print(f"    {k}={env[k]}")
            continue

        env       = build_env(exp, alpha, defaults, data_paths)
        job_start = time.time()

        eval_status = None
        try:
            result = subprocess.run(
                [sys.executable, "train_single.py"],
                env=env,
                check=False,
            )
            elapsed = time.time() - job_start
            status  = "OK" if result.returncode == 0 else f"FAILED (rc={result.returncode})"
            if result.returncode != 0:
                failed.append(label)
            elif not args.no_eval_test and exp.get("eval_test", defaults.get("eval_test", True)):
                print(f"  Evaluating on test split...")
                eval_status = eval_test_set(exp, alpha, ldir, env, defaults)
                print(f"  {eval_status}")
        except KeyboardInterrupt:
            print("\n[Interrupted — stopping runner]")
            break
        except Exception as exc:
            elapsed = time.time() - job_start
            status  = f"ERROR: {exc}"
            failed.append(label)

        line = (
            f"{datetime.now().isoformat()}  {status}"
            f"  elapsed={hms(elapsed)}"
            + (f"  {eval_status}" if eval_status else "")
            + f"  {label}"
        )
        print(f"\n  {line}")
        with open(log_path, "a") as lf:
            lf.write(line + "\n")

    # ── Summary ───────────────────────────────────────────────────────────────
    total_elapsed = time.time() - global_start
    print(f"\n{'━'*65}")
    print(f"  All jobs finished in {hms(total_elapsed)}")
    if failed:
        print(f"  FAILED ({len(failed)}):")
        for f_ in failed:
            print(f"    • {f_}")
    else:
        print("  All jobs completed successfully ✓")
    print(f"{'━'*65}\n")

    if not args.dry_run:
        with open(log_path, "a") as lf:
            lf.write(
                f"Run finished : {datetime.now().isoformat()}"
                f"  total={hms(total_elapsed)}\n"
            )


if __name__ == "__main__":
    main()

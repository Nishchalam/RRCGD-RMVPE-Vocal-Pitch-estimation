"""
train_single.py  (ACS of RRCGD)
────────────────────────────────
Runs ONE experiment (one alpha value, one dataset/fmax/window_length/
bin_method combination).  All parameters come from environment variables set
by run_experiments.py before the subprocess is spawned.

Supported datasets
──────────────────
  mir1k  → MIR1K  class, SR = 16 000 Hz  (music, .pv pitch labels)
  ptdb   → PTDB   class, SR =  8 000 Hz  (speech, .npy pitch labels)

SR flows through the env-driven constant so both dataset classes and the
model constructor receive the correct value automatically:

    EXP_SAMPLE_RATE=8000  →  constants.SAMPLE_RATE=8000
      → PTDB.HOP_LENGTH  = int(20/1000 * 8000)  = 160 samples
      → PTDB.seq_len     = int(2.55   * 8000)   = 20400 samples
      → E2E(hop_samples) = E2E(160, ...)

Direct invocation example (PTDB):
    EXP_SAMPLE_RATE=8000 EXP_F_MAX=4000 EXP_WINDOW_LENGTH=2048 \\
    EXP_BIN_METHOD=linear_fmax                                   \\
    EXP_DATASET=ptdb EXP_DATATYPE=full_data                      \\
    EXP_DATA_PATH="/path/to/PTDB_data"                           \\
    EXP_TRAIN_SPLITS="['train']" EXP_VALID_SPLITS="['valid']"    \\
    EXP_ALPHA=16 EXP_GAMMA=0                                     \\
    python train_single.py
"""

import os
import ast

# ── All parameters from env vars ─────────────────────────────────────────────
# EXP_SAMPLE_RATE must be read before importing src so constants.py sees it.
DATASET      = os.environ["EXP_DATASET"]
DATATYPE     = os.environ["EXP_DATATYPE"]
DATA_PATH    = os.environ["EXP_DATA_PATH"]
TRAIN_SPLITS = ast.literal_eval(os.environ["EXP_TRAIN_SPLITS"])
VALID_SPLITS = ast.literal_eval(os.environ["EXP_VALID_SPLITS"])
ALPHA        = float(os.environ["EXP_ALPHA"])
GAMMA        = float(os.environ.get("EXP_GAMMA",          "0"))
LR           = float(os.environ.get("EXP_LR",             "1e-3"))
BATCH_SIZE   = int(os.environ.get("EXP_BATCH_SIZE",        "8"))
SEQ_L        = float(os.environ.get("EXP_SEQ_L",           "2.55"))
HOP_LENGTH   = int(os.environ.get("EXP_HOP_LENGTH",        "20"))
CLIP_GRAD    = int(os.environ.get("EXP_CLIP_GRAD_NORM",    "3"))
CLEAN_ONLY   = os.environ.get("EXP_CLEAN_ONLY", "0") == "1"
USE_CACHE    = os.environ.get("EXP_USE_CACHE", "1") == "1"
SEED          = int(os.environ.get("EXP_SEED",           "42"))
DETERMINISTIC = os.environ.get("EXP_DETERMINISTIC", "1") == "1"
# Opt-in escape hatch for PreRepro-replication runs: skips set_seed() and
# the DataLoader's manual generator/worker_init_fn entirely, so RNG state
# is whatever plain unseeded PyTorch/NumPy/random give you -- matching
# the PreRepro codebase (rrcgd_rmvpe_final.zip), which never called
# torch.manual_seed or any seeding helper. Default stays "0" so every
# existing invocation (the live pipeline included) is unaffected.
NO_SEED       = os.environ.get("EXP_NO_SEED", "0") == "1"
PATIENCE_EPOCHS = int(os.environ.get("EXP_PATIENCE_EPOCHS", "4"))
# See run_experiments.py's --repeat help text: distinguishes multiple runs
# of the identical (seed, config) combo under deterministic=False, where
# they'd otherwise collide on the same run dir. Default 1 (untagged) so
# every existing invocation is unaffected.
REPEAT        = int(os.environ.get("EXP_REPEAT", "1"))

# ─────────────────────────────────────────────────────────────────────────────
# src imports come AFTER env vars are set so constants.py reads EXP_SAMPLE_RATE
import torch
from torch import nn
from torch.nn.utils import clip_grad_norm_
from torch.optim.lr_scheduler import StepLR
from torch.utils.data import DataLoader
from torch.utils.tensorboard import SummaryWriter
from torch.cuda.amp import autocast, GradScaler
from tqdm import tqdm
import numpy as np

from src import MIR1K, PTDB, E2E, cycle, SAMPLE_RATE, FL, F_MAX, WINDOW_LENGTH, PITCH_TH,N_BINS, set_seed, seed_worker
from src.constants import BIN_METHOD
from evaluate import evaluate, evaluate_batch


# ── Dataset class registry ────────────────────────────────────────────────────
_DATASET_CLS = {
    "mir1k": MIR1K,
    "ptdb":  PTDB,
}

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
}


def train():
    # Before anything that draws from an RNG — model init and DataLoader
    # shuffling both depend on this having already run.
    if not NO_SEED:
        set_seed(SEED, deterministic=DETERMINISTIC)

    device  = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # fp32 by default: fp16 autocast made high-alpha runs collapse to non-finite loss
    # (no checkpoint); EXP_AMP=1 restores the old fp16 behaviour.
    use_amp = device.type == "cuda" and os.environ.get("EXP_AMP", "0") == "1"

    # fp32 jobs need ~50% more GPU memory than fp16 ones; wait for headroom before
    # allocating so concurrent jobs can't push each other into OOM.
    if device.type == "cuda":
        import subprocess as _sp, time as _t
        _gpu = os.environ.get("CUDA_VISIBLE_DEVICES", "0").split(",")[0]
        _need = (9500 if not use_amp else 6500) if int(os.environ.get("EXP_N_BINS", 256)) >= 512 \
                else (6500 if not use_amp else 4500)
        while True:
            try:
                _o = _sp.run(["nvidia-smi", "-i", _gpu, "--query-gpu=memory.used,memory.total",
                              "--format=csv,noheader,nounits"], capture_output=True, text=True).stdout
                _u, _tot = (int(x) for x in _o.split(","))
            except Exception:
                break
            if _tot - _u >= _need:
                break
            print(f"[wait] GPU{_gpu} free {_tot - _u} MiB < {_need} MiB needed; waiting", flush=True)
            _t.sleep(60)

    # ── Validate dataset choice ───────────────────────────────────────────────
    if DATASET not in _DATASET_CLS:
        raise ValueError(
            f"Unknown EXP_DATASET='{DATASET}'. "
            f"Valid options: {list(_DATASET_CLS.keys())}"
        )
    DatasetCls = _DATASET_CLS[DATASET]

    # ── Log / checkpoint dir ──────────────────────────────────────────────────
    # Include SR in the tag so ptdb and mir1k log dirs never collide even if
    # they share the same f_max / window_length / bin_method combination.
    dataset_tag = f"{DATASET}_sr{SAMPLE_RATE}"

    bins_tag = f"bins{N_BINS}"

    fmax_tag = (
        f"{F_MAX}_FL_{WINDOW_LENGTH}_{bins_tag}"
    )

    method_tag = _METHOD_TAGS.get(
        BIN_METHOD,
        BIN_METHOD
    )

    # Only tag non-default seeds so every existing run path (all seeded at
    # the project default, 42) stays exactly as-is -- same backward-compat
    # convention already used for window-length in export_validation_csv.py.
    # Without this, two runs that differ ONLY by seed would silently share
    # one directory: --skip-existing would treat the second as already
    # done, or worse, a second training pass could overwrite the first
    # seed's checkpoint mid-run.
    seed_tag = "" if SEED == 42 else f"_seed{SEED}"
    # no_seed wins over seed_tag -- see run_experiments.py's matching
    # comment; must stay byte-for-byte in sync with logdir_for() there.
    if NO_SEED:
        seed_tag = "_noseed"
    # See run_experiments.py's logdir_for() matching comment -- must stay
    # byte-for-byte in sync with it.
    repeat_tag = "" if REPEAT == 1 else f"_rep{REPEAT}"
    runs_root = os.environ.get("EXP_RUNS_ROOT") or ("runs_seeded" if DETERMINISTIC else "runs_deterministic_false")

    parent_dir = (
        f"{runs_root}/pitch_th_{PITCH_TH}"
        f"/{dataset_tag}"
        f"/{DATATYPE}"
        f"/{fmax_tag}"
        f"/{method_tag}_{LR}_{BATCH_SIZE}{seed_tag}{repeat_tag}"
    )

    loss_tag = f"Pitch_FL_{ALPHA}_{GAMMA}"

    logdir = f"{parent_dir}/{loss_tag}"

    os.makedirs(logdir, exist_ok=True)
    writer = SummaryWriter(logdir)

    # Hop length in samples — correct for this dataset's SR
    hop_samples = int(HOP_LENGTH / 1000 * SAMPLE_RATE)

    print(f"\n{'='*65}")
    print(f"  Experiment  : {logdir}")
    print(f"  Dataset cls : {DatasetCls.__name__}  (SR={SAMPLE_RATE} Hz)")
    print(f"  Split       : {DATATYPE}")
    print(f"  F_MAX       : {F_MAX} Hz   WL : {WINDOW_LENGTH}")
    print(f"  Bin method  : {BIN_METHOD} → '{method_tag}'")
    print(f"  Hop         : {HOP_LENGTH} ms = {hop_samples} samples")
    print(f"  Alpha       : {ALPHA}   Gamma : {GAMMA}")
    print(f"  Data cond.  : {'clean only' if CLEAN_ONLY else 'clean + noisy'} | cache={USE_CACHE}")
    print(f"  LR          : {LR}   BS : {BATCH_SIZE}")
    print(f"  Device      : {device}   AMP : {use_amp}")
    if NO_SEED:
        print(f"  Seed        : NONE (unseeded, PreRepro-replication mode)")
    else:
        print(f"  Seed        : {SEED}   Deterministic : {DETERMINISTIC}")
    print(f"{'='*65}\n")

    # ── Datasets ──────────────────────────────────────────────────────────────
    # PTDB and MIR1K share the same (path, hop_length, sequence_length, groups,
    # use_cache) constructor signature and caching design.
    train_dataset = DatasetCls(
        DATA_PATH, HOP_LENGTH, SEQ_L, TRAIN_SPLITS, use_cache=USE_CACHE
    )
    validation_dataset = DatasetCls(
        DATA_PATH, HOP_LENGTH, None, VALID_SPLITS, use_cache=USE_CACHE
    )

    print(f"Train   : {len(train_dataset)} samples")
    print(f"Valid   : {len(validation_dataset)} samples")

    if NO_SEED:
        # No generator, no worker_init_fn -- shuffle order and each
        # DataLoader worker's RNG are whatever plain PyTorch gives you,
        # matching PreRepro's DataLoader construction exactly.
        data_loader = DataLoader(
            train_dataset, BATCH_SIZE,
            shuffle=True, drop_last=True,
            num_workers=4, pin_memory=True, persistent_workers=True,
        )
    else:
        # Generator fixes the shuffle order; seed_worker fixes each of the 4
        # worker subprocesses' own RNGs, which torch does not derive from the
        # global seed on its own.
        g = torch.Generator()
        g.manual_seed(SEED)

        data_loader = DataLoader(
            train_dataset, BATCH_SIZE,
            shuffle=True, drop_last=True,
            num_workers=4, pin_memory=True, persistent_workers=True,
            worker_init_fn=seed_worker, generator=g,
        )
    if (DATASET == "mir1k"):
        validation_interval       = len(data_loader)
        iterations                = len(data_loader) * 50
        learning_rate_decay_steps = len(data_loader) * 5
    else:
        validation_interval       = len(data_loader)
        iterations                = len(data_loader) * 15
        learning_rate_decay_steps = len(data_loader) * 5
    # ── Model / optimiser / scaler ────────────────────────────────────────────
    model = nn.DataParallel(
        E2E(hop_samples, 4, 1, (2, 2), device=device)
    ).to(device)
    print(f"  Input bins  : {N_BINS}")
    print(f"  Effective bins : {model.module.rrcgd.n_bins}")
    print(f"  Raw retained bins : {model.module.rrcgd.raw_n_bins}")
    print(f"  Frequency padding : {model.module.rrcgd.pad_amount}")
    optimizer = torch.optim.Adam(model.parameters(), LR)
    scheduler = StepLR(optimizer, step_size=learning_rate_decay_steps, gamma=0.98)
    scaler    = GradScaler(enabled=use_amp)

    # ── Training loop ─────────────────────────────────────────────────────────
    loop = tqdm(range(1, iterations + 1))
    loop.set_description(f"{method_tag} {fmax_tag} | {loss_tag}")

    RPA, RCA, OA, VFA, VR, it = 0, 0, 0, 0, 0, 0

    for i, data in zip(loop, cycle(data_loader)):
        roots       = data["roots"].to(device, non_blocking=True)
        pitch_label = data["pitch"].to(device, non_blocking=True)

        with autocast(enabled=use_amp):
            pitch_pred = model(roots)
            T          = pitch_label.shape[1]
            pitch_pred = pitch_pred[:, :T, :]
            loss       = FL(pitch_pred, pitch_label, ALPHA, GAMMA)

        if not torch.isfinite(loss):
            # RRCGD is now clamped so this shouldn't trigger, but a skipped
            # step here is cheap insurance against another silent BatchNorm
            # corruption (GradScaler only guards the gradient path, not the
            # BN running-stat update baked into this forward pass).
            optimizer.zero_grad(set_to_none=True)
            # A near-resonance episode forces _chunked_rrcgd_sum's chunk
            # tensors into shapes that don't match the caching allocator's
            # existing free blocks, so its reserved (not just allocated)
            # pool grows with every skipped iteration and never shrinks on
            # its own -- across thousands of consecutive skips (the common
            # case once one of these episodes starts) that reserved pool is
            # what pushes nvidia-smi's "used" figure up, not a real leak.
            # Periodically returning it lets a concurrent job on the same
            # GPU keep its headroom instead of this job slowly starving it.
            del loss, pitch_pred
            if i % 50 == 0:
                torch.cuda.empty_cache()
            print(f"\n[warn] non-finite loss at iter {i}, skipping step")
            continue

        optimizer.zero_grad(set_to_none=True)
        scaler.scale(loss).backward()

        if CLIP_GRAD:
            scaler.unscale_(optimizer)
            clip_grad_norm_(model.parameters(), CLIP_GRAD)

        scaler.step(optimizer)
        scaler.update()
        scheduler.step()

        writer.add_scalar("loss/loss_pitch", loss.item(), global_step=i)

        with torch.no_grad():
            train_metrics = evaluate_batch(pitch_pred, pitch_label, HOP_LENGTH)
        writer.add_scalar("train/RPA", train_metrics["RPA"], global_step=i)
        writer.add_scalar("train/OA",  train_metrics["OA"],  global_step=i)
        writer.add_scalar("train/VR",  train_metrics["VR"],  global_step=i)

        loop.set_postfix(loss=f"{loss.item():.4f}", best_rpa=f"{RPA:.4f}")

        # ── Validation ────────────────────────────────────────────────────────
        if i % validation_interval == 0:
            model.eval()
            with torch.no_grad():
                metrics = evaluate(
                    validation_dataset, model.module, HOP_LENGTH, device
                )
                for key, value in metrics.items():
                    writer.add_scalar(
                        "stage_pitch/" + key, np.mean(value), global_step=i
                    )

                rpa = np.mean(metrics["RPA"])
                rca = np.mean(metrics["RCA"])
                oa  = np.mean(metrics["OA"])
                vr  = np.mean(metrics["VR"])
                vfa = np.mean(metrics["VFA"])

                print(
                    f"\nIter {i}: RPA={rpa:.4f}, RCA={rca:.4f}, "
                    f"OA={oa:.4f}, VR={vr:.4f}, VFA={vfa:.4f}"
                )

                if rpa >= RPA:
                    RPA, RCA, OA, VR, VFA, it = rpa, rca, oa, vr, vfa, i

                    with open(os.path.join(logdir, "result.txt"), "a") as f:
                        f.write(
                            f"{i}\t{RPA:.4f}\t{RCA:.4f}\t"
                            f"{OA:.4f}\t{VR:.4f}\t{VFA:.4f}\n"
                        )

                    torch.save(
                        model.module,
                        os.path.join(logdir, f"model-1-{i}.pt"),
                    )
                    # torch.save(
                    #     optimizer.state_dict(),
                    #     os.path.join(logdir, "last-optimizer-state.pt"),
                    # )
                    print("✓ New best model saved")

                    ckpts = sorted(
                        [f for f in os.listdir(logdir)
                         if f.startswith("model-1-") and f.endswith(".pt")],
                        key=lambda x: int(x.split("-")[-1].split(".")[0]),
                    )
                    for ckpt in ckpts[:-1]:
                        os.remove(os.path.join(logdir, ckpt))

            model.train()

        if i - it > len(data_loader) * PATIENCE_EPOCHS:
            # Best checkpoints in this sweep typically land by epoch 3-14/30
            # (clean) or 1-4/30 (clean_noisy); patience defaults to 4 epochs
            # (was 8, which meant ~70%+ of a run's wall-clock was spent
            # idling after the actual best checkpoint -- costly at low batch
            # sizes where each epoch is len(data_loader) steps) but is
            # configurable via EXP_PATIENCE_EPOCHS for tracks that need more
            # room, e.g. the 128ms window track using 8 to rule out
            # premature stopping as a factor in its lower-than-expected RPA.
            print(f"\nEarly stopping triggered (patience={PATIENCE_EPOCHS} epochs).")
            break

    writer.close()
    print(f"\nFinished: {logdir}  |  best RPA = {RPA:.4f} @ iter {it}\n")


if __name__ == "__main__":
    train()

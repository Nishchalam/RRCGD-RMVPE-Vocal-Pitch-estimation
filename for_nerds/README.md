# For nerds: auto GPU scheduler

A single script that runs an entire sweep config unattended across
however many GPUs the machine actually has, without you hand-picking
which experiment goes on which GPU.

```bash
./auto_gpu_scheduler.sh ../configs/experiments_clean_256bins_bs8.yaml
./auto_gpu_scheduler.sh ../configs/experiments_clean_noisy_512bins_bs8.yaml 9500
```

## What it does

1. **Detects GPUs and capacity.** Reads GPU count and per-GPU total
   memory from `nvidia-smi`, and divides by `per_job_mib` (second
   argument, default 6500) to decide how many jobs each GPU can hold
   at once.
2. **Launches one screen per experiment**, in the order the config
   lists them, onto whichever GPU currently has a free slot *and*
   enough actually-free memory (re-checked immediately before every
   launch — never batch-launches off a stale reading). Each screen runs
   `run_experiments.py --only <experiment> --skip-existing`, which
   itself sequentially sweeps that experiment's full alpha list.
3. **Resumable at two levels.** `--skip-existing` means a
   killed/restarted scheduler run won't re-launch an experiment whose
   full alpha sweep already finished; and as of this version,
   `train_single.py` itself auto-resumes mid-experiment too (model +
   optimizer + scheduler + scaler state are all saved alongside the
   best checkpoint — see the main [INSTALL.md](../INSTALL.md)) — so a
   job killed partway through one alpha picks back up from its last
   checkpoint instead of restarting that alpha from iteration 0.
4. **Auto-exports results every 10 completed jobs** (counted by
   on-disk `test_eval.csv` files, never by parsing logs — log-based
   progress tracking is a known footgun here, see the comments in
   `export_validation_csv.py`): re-runs `export_validation_csv.py` and
   `compile_best_results.py` so the aggregated CSVs stay current
   through a long unattended run.
5. **Exits cleanly** once every experiment in the config has a
   finished screen, after one final CSV export.

## Picking `per_job_mib`

The project's own empirically-measured per-job memory (fp32, this
project's batch size):

| Bins | Technique | ~Memory |
|---|---|---|
| 256 | full-range (E1/E2) | ~6.5GB |
| 256 | frequency/bin-truncated (E3/E4) | ~3.8GB |
| 512 | full-range (E1/E2) | ~9.5GB |
| 512 | frequency/bin-truncated (E3/E4) | ~4.5GB |

The default (6500) is a reasonable single number for a 256-bin sweep
with a mix of techniques; pass a larger value for a 512-bin config, or
a smaller one if you know your config is all truncated techniques and
want more concurrency. Err high rather than low — this is a capacity
*and* a launch-gate number, so underestimating risks OOM, while
overestimating just under-uses the card.

## Logs

Per-experiment stdout goes to `logs/<experiment>.log` (gitignored, not
committed) next to this script. The scheduler's own progress prints to
stdout — run it inside `screen`/`tmux` for a long sweep:

```bash
screen -dmS scheduler bash -c './auto_gpu_scheduler.sh ../configs/experiments_clean_256bins_bs8.yaml'
screen -r scheduler   # attach to watch; Ctrl+A D to detach
```

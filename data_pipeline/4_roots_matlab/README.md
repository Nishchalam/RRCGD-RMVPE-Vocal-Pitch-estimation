# Step 4: Polynomial root extraction (MATLAB)

A patched, Linux/Mac-portable build of the **Lindsey-Fox (WebLF21)**
fast polynomial root finder, used here to extract the per-frame
polynomial roots that `src/spec.py`'s `RRCGDSpectrogram` consumes.

The original WebLF21 package ships several bugs and Windows-only
assumptions that fail out of the box on modern MATLAB / non-Windows
systems. What's here is the patched version — see
[FIXES_SUMMARY.md](FIXES_SUMMARY.md) for the full list of what was
broken and how it was fixed.

## Third-party code notice

This directory bundles the Lindsey-Fox toolbox and Zhonggang Zeng's
`multroot` package. `multroot`'s own header states it is "freely
released for research exchange only" — verify the licensing terms of
each bundled component yourself before any non-research use. This code
is included for reproducibility of the root-extraction step; it is not
original work of this project, and **it is explicitly excluded from
the repository's MIT [LICENSE](../../LICENSE)** — it retains whatever
terms its own authors gave it.

## Layout

| Path | Purpose |
|---|---|
| `lroots/`, `shared/`, `Unwrap/`, `util/`, `broots/`, `multroot/`, `simplified/` | The Lindsey-Fox toolbox itself (patched) |
| `mexFiles/` | Optional compiled accelerators (recompile with `mex *.c` for a real speedup on Linux/Mac) |
| `setup_lf_path.m` | Run once per MATLAB session before calling `lroots`/`broots`/`unwrapz` |
| `extract_roots_batch.m` | Main batch driver: frames a directory of `.wav` files, windows each frame, computes polynomial roots per frame |
| `run_experiment.m` | Calls `extract_roots_batch` over PTDB's train/valid/test splits at 32ms/10ms (this project's setting) |
| `monitor_and_halt.sh` | Watches system load average and kills the MATLAB root-extraction job if load stays too high for too long |

## Requirements

- MATLAB with the Parallel Computing Toolbox (`extract_roots_batch.m` uses `parfor`/`parcluster`)
- Audio files readable by `audioread` (`.wav`) — i.e. the output of
  [step 2](../2_split/) and [step 3](../3_augment/)

## Usage

### 1. Set up the path (once per MATLAB session)

```matlab
setup_lf_path
rootsIni
```

### 2. Run the project's extraction set

```matlab
setup_lf_path
rootsIni
run_experiment
```

This extracts 32ms-window, 10ms-hop roots for `train/`, `valid/`, and
`test/`, writing one `.mat` file per input `.wav` into an auto-named
`roots_WL_32ms_HOP_10ms/` subfolder next to the audio. Already-processed
files are skipped automatically, so a killed/restarted run only
processes what's left.

For a long batch run, launch MATLAB headless inside a `screen` session
so it survives disconnects:

```bash
screen -dmS roots_extraction matlab -nodisplay -nosplash -nodesktop -r "setup_lf_path; rootsIni; run_experiment; exit"
screen -r roots_extraction   # attach to watch progress; Ctrl+A D to detach
```

### 3. (Optional) Extract roots from any other directory directly

```matlab
% extract_roots_batch(audio_dir, sample_rate_hz, frame_size_ms, hop_length_ms, method)
extract_roots_batch('/path/to/wavs', 8000, 32, 10, 'lroots');
```

`method` is `'lroots'` (Lindsey-Fox, default) or `'roots'` (MATLAB's
built-in, for comparison).

### 4. Guard against overload while it runs

```bash
./monitor_and_halt.sh [LOAD_THRESHOLD] [POLL_SECONDS] [SUSTAIN_CHECKS]
```

Defaults: threshold `20`, poll every `30s`, requires `3` consecutive
over-threshold checks before it sends `SIGTERM`/`SIGKILL` to the MATLAB
root-extraction job specifically (matched precisely so it won't touch
unrelated jobs). It exits on its own once the job finishes.

## Notes

- Polynomial degree scales with window length at a fixed sample rate.
  At 8kHz/32ms (this project's setting) the degree is small enough that
  the pure-MATLAB `lroots` path is fine without the compiled mex
  accelerators; see [FIXES_SUMMARY.md](FIXES_SUMMARY.md) for
  measured/extrapolated timing at larger window lengths.
- Output goes next to the input audio (`roots_WL_<ms>ms_HOP_<ms>ms/`),
  which is exactly where `src/dataset.py` expects to find it — no
  further file moving needed before step 6.

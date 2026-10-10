# Step 1: Audio download and environment setup

## 1a. Download PTDB-TUG

```bash
./download_ptdb.sh /path/to/PTDB_raw
```

This fetches the single ~3.9GB zip from TU Graz's SPSC database page
(no registration required), verifies its MD5 checksum, and extracts
it. If the script's URL ever goes stale, get the file manually from:
https://www2.spsc.tugraz.at/databases/PTDB-TUG/

Expected layout after extraction:

```
PTDB_raw/
  SPEECH DATA/
    FEMALE/
      F01/ ... F10/
    MALE/
      M01/ ... M10/
  DOCUMENTATION/
```

Each speaker folder contains paired microphone recordings
(`mic_<speaker>_<utterance>.wav`) and laryngograph-derived pitch
reference files (`ref_<speaker>_<utterance>.f0`), though the exact
subfolder nesting under each speaker can vary by release — [step
2](../2_split/split_ptdb.py)'s script searches recursively for these
rather than assuming a fixed depth.

## 1b. Set up the Python/MATLAB environment

Everything else (conda env, PyTorch install, Python package
requirements) is covered in the top-level [INSTALL.md](../../INSTALL.md)
— the environment name used throughout this project is `gd_exps`. You
also need MATLAB with the Parallel Computing Toolbox for [step
4](../4_roots_matlab/) (root extraction); there is no way around that
requirement for this feature.

```bash
conda env create -f ../../environment.yml -n gd_exps
conda activate gd_exps
pip install torch==2.10.0 torchaudio==2.10.0 torchvision==0.25.0 \
    --index-url https://download.pytorch.org/whl/cu128
pip install -r ../../requirements.txt
```

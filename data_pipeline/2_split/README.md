# Step 2: Speaker-disjoint train/valid/test split

```bash
python split_ptdb.py --raw-dir /path/to/PTDB_raw --out dataset/PTDB_data_10ms_hop --dry-run
python split_ptdb.py --raw-dir /path/to/PTDB_raw --out dataset/PTDB_data_10ms_hop
```

PTDB-TUG has 20 speakers (F01–F10, M01–M10). This project's split is
**speaker-disjoint** — no speaker's audio appears in more than one of
train/valid/test, which is the whole point of this project's headline
result (the RRCGD paper's own 96% number was measured under a *leaky*
split with shared speakers across train/test; under this disjoint
split the honest ceiling is ~90–94%).

| Split | Speakers | Count |
|---|---|---|
| train | F01–F06, M01–M06 | 12 |
| valid | F07, F08, M07, M08 | 4 |
| test | F09, F10, M09, M10 | 4 |

Edit `SPEAKER_SPLIT` in `split_ptdb.py` if you want a different split
— just keep it speaker-disjoint.

Output: `.wav` copied as-is (no resampling here — [step
4](../4_roots_matlab/)'s MATLAB extractor resamples to the project's
8kHz internally), plus a `.npy` per utterance holding the `.f0`
reference's pitch column (Hz) as a 1D array, frame-aligned with the
10ms-hop roots step 4 will compute from the same audio.

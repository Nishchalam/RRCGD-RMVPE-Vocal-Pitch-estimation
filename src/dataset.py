'''
This is dataset.py from where RMPVE takes rrcgd features and trains on the proper data split criterion made. Clause to be added for taking in just clean or clean+noisy while training in here
'''
import os
import re
import librosa
import numpy as np
import torch
import scipy.io #for importing .mat root files
from torch.utils.data import Dataset
from tqdm import tqdm
from glob import glob
import pandas as pd
from .constants import *


def _atomic_torch_save(obj, path):
    """torch.save that's safe when two training processes (e.g. a mel and
    a linear job sharing the same window-length/hop cache dir, since the
    cache key doesn't depend on bin_method) race to populate the same
    cold cache path: each writer saves to its own pid+object-id-tagged
    temp file, then os.replace() atomically installs it. os.replace is
    a single rename syscall on the same filesystem, so a reader can never
    observe a half-written file -- unlike torch.save(obj, path) directly,
    which two concurrent writers can interleave into a truncated zip that
    fails with "failed finding central directory" forever after (the
    corrupt file's mere presence short-circuits future cache-hit checks).
    """
    tmp_path = f"{path}.tmp.{os.getpid()}.{id(obj)}"
    torch.save(obj, tmp_path)
    os.replace(tmp_path, path)


def _atomic_write_lines(path, lines):
    """Same rationale as _atomic_torch_save, for the plain-text manifest
    files written alongside chunked caches."""
    tmp_path = f"{path}.tmp.{os.getpid()}"
    with open(tmp_path, 'w') as f:
        for line in lines:
            f.write(line + '\n')
    os.replace(tmp_path, path)


class MIR1K(Dataset):
    """
    MIR-1K dataset for the optimized RRCGD-RMVPE experiment (singer-disjoint
    new_data_split/, polyphonic 'audio' + isolated 'clean_vocals' groups,
    each with clean + 4 SNR-augmented noise levels).

    Mirrors PTDB's design so both datasets share one mental model:
      1. Uses roots_WL_<window>ms_HOP_<hop>ms, same naming convention as
         PTDB's root-extraction directories (e.g. roots_WL_32ms_HOP_20ms).
      2. clean-only vs clean+noisy is controlled by EXP_CLEAN_ONLY.
      3. Noisy roots (snrNN_<stem>.mat, label snrNN_<stem>.f0) reuse the
         clean <stem>.pv label by stripping the SNR prefix -- verified
         byte-identical to the corresponding snrNN_<stem>.f0 on disk.
      4. .mat parsing is cached to .pt files, exactly like PTDB.
      5. condition_at() exposes clean/snrNN per item for evaluate_mir1k_checkpoint.py.

    Label format differs from PTDB: MIR-1K's .pv/.f0 files store MIDI note
    numbers (0 = unvoiced), not Hz, so _parse() converts via the standard
    MIDI-to-Hz formula before going to cents (PTDB's .npy labels are already Hz).

    Frame-length convention differs by group: clean-vocal-stem roots have one
    fewer invalid boundary frame than polyphonic-mix roots for the same window
    length (WINDOW_LENGTH-1 vs WINDOW_LENGTH-2) -- preserved from the original
    loader, which established this by inspecting the actual .mat frame counts.
    """

    def __init__(self, path, hop_length=None, sequence_length=None, groups=None,
                 use_cache=None):
        self.path = path
        hop_ms = HOP_MS if hop_length is None else int(hop_length)
        self.HOP_LENGTH = int(round(hop_ms / 1000 * SAMPLE_RATE))
        self.hop_ms = hop_ms
        self.seq_len = None if not sequence_length else int(sequence_length * SAMPLE_RATE)
        self.num_class = N_CLASS
        self.groups = list(groups or [])
        self.use_cache = USE_CACHE if use_cache is None else bool(use_cache)
        self.chunked = self.seq_len is not None

        if self.chunked:
            self.n_steps = self.seq_len // self.HOP_LENGTH + 1
        else:
            self.n_steps = None

        print(
            f"Indexing {len(self.groups)} group{'s' if len(self.groups) != 1 else ''} "
            f"of MIR1K at {path} | window={WINDOW_LENGTH} samples "
            f"({WINDOW_LENGTH / SAMPLE_RATE * 1000:.0f} ms), hop={self.hop_ms} ms "
            f"({self.HOP_LENGTH} samples), seq_frames={self.n_steps}, "
            f"CLEAN_ONLY={CLEAN_ONLY}, cache={'on' if self.use_cache else 'off'}"
        )

        self.index = []
        self.root_paths = []
        self.condition = []
        skipped = 0

        for group in self.groups:
            pairs = self.files(group)
            for root_path, label_path in tqdm(
                pairs, desc=f"Indexing group {group}"
            ):
                if self.use_cache:
                    cache_entries = self._ensure_cache(group, root_path, label_path)
                    if not cache_entries:
                        skipped += 1
                        continue
                    for cp in cache_entries:
                        self.index.append(cp)
                        self.root_paths.append(root_path)
                        self.condition.append(parse_condition(root_path))
                else:
                    self.index.append((root_path, label_path))
                    self.root_paths.append(root_path)
                    self.condition.append(parse_condition(root_path))

        if skipped:
            print(f"[INFO] Skipped {skipped} file(s) with no usable frames.")
        print(f"Indexed {len(self.index)} dataset item(s) for {self.groups}.")

    @staticmethod
    def availabe_groups():
        return ['train_audio', 'valid_audio', 'test_audio',
                 'train_clean_vocals', 'valid_clean_vocals', 'test_clean_vocals']

    def files(self, group):
        window_ms = int(round(WINDOW_LENGTH / SAMPLE_RATE * 1000))
        roots_dir = os.path.join(
            self.path, group,
            f"roots_WL_{window_ms}ms_HOP_{self.hop_ms}ms"
        )
        if not os.path.isdir(roots_dir):
            raise FileNotFoundError(
                f"Expected MIR1K roots directory not found:\n  {roots_dir}\n"
                f"Check SAMPLE_RATE={SAMPLE_RATE}, WINDOW_LENGTH={WINDOW_LENGTH}, "
                f"HOP_MS={self.hop_ms}."
            )

        root_files = sorted(glob(os.path.join(roots_dir, '*.mat')))
        if CLEAN_ONLY:
            root_files = [
                rf for rf in root_files
                if strip_snr_prefix(os.path.splitext(os.path.basename(rf))[0])[1] is None
            ]

        pairs = []
        missing = []
        for root_path in root_files:
            stem = os.path.splitext(os.path.basename(root_path))[0]
            clean_stem, _ = strip_snr_prefix(stem)
            label_path = os.path.join(self.path, group, clean_stem + '.pv')
            if not os.path.isfile(label_path):
                missing.append((stem, label_path))
            else:
                pairs.append((root_path, label_path))

        if missing:
            preview = "\n".join(f"  {s} -> {p}" for s, p in missing[:5])
            raise FileNotFoundError(
                f"{len(missing)} MIR1K root file(s) have no matching clean label.\n{preview}"
            )
        return pairs

    def _cache_dir(self, group):
        window_ms = int(round(WINDOW_LENGTH / SAMPLE_RATE * 1000))
        return os.path.join(
            self.path, group,
            f"_frame_cache_WL{window_ms}ms_HOP{self.hop_ms}ms_FREF10.0_NCLASS{N_CLASS}"
        )

    def _cache_stem(self, root_path):
        return os.path.splitext(os.path.basename(root_path))[0]

    def _ensure_cache(self, group, root_path, label_path):
        """Create cache once and return a list of cache files for this item."""
        cache_dir = self._cache_dir(group)
        os.makedirs(cache_dir, exist_ok=True)
        stem = self._cache_stem(root_path)

        if self.chunked:
            manifest = os.path.join(cache_dir, stem + ".manifest")
            if os.path.isfile(manifest):
                entries = []
                with open(manifest) as f:
                    for line in f:
                        cp = line.strip()
                        if cp and os.path.isfile(cp):
                            entries.append(cp)
                if entries:
                    return entries

            item = self._parse(group, root_path, label_path)
            if item is None:
                return []

            T = item['roots'].shape[0]
            entries = []
            n_full = T // self.n_steps
            for i in range(n_full):
                begin = i * self.n_steps
                end = begin + self.n_steps
                cp = os.path.join(cache_dir, f"{stem}__chunk{i:04d}.pt")
                if not os.path.isfile(cp):
                    _atomic_torch_save({
                        'roots': item['roots'][begin:end],
                        'pitch': item['pitch'][begin:end],
                        'voice': item['voice'][begin:end],
                        'file': root_path,
                    }, cp)
                entries.append(cp)

            if T % self.n_steps != 0:
                cp = os.path.join(cache_dir, f"{stem}__chunk{n_full:04d}.pt")
                if not os.path.isfile(cp):
                    _atomic_torch_save({
                        'roots': item['roots'][-self.n_steps:],
                        'pitch': item['pitch'][-self.n_steps:],
                        'voice': item['voice'][-self.n_steps:],
                        'file': root_path,
                    }, cp)
                entries.append(cp)

            _atomic_write_lines(manifest, entries)
            del item
            return entries

        # Full-file cache for validation/test.
        cp = os.path.join(cache_dir, stem + ".pt")
        if os.path.isfile(cp):
            return [cp]
        item = self._parse(group, root_path, label_path)
        if item is None:
            return []
        _atomic_torch_save(item, cp)
        del item
        return [cp]

    def _parse(self, group, root_path, label_path):
        """Parse one .mat + .pv/.f0 pair into CPU tensors."""
        mat_data = scipy.io.loadmat(root_path)
        if 'roots' not in mat_data:
            raise KeyError(f"{root_path} does not contain MATLAB variable 'roots'")
        raw_roots = mat_data['roots']

        # Verified against actual roots_WL_*ms_HOP_*ms .mat files: both
        # 'audio' (polyphonic) and 'clean_vocals' (isolated vocal stem)
        # groups produce WINDOW_LENGTH-1 roots per frame under this root-
        # extraction scheme -- same convention as PTDB. (The old commented-
        # out loader's audio/clean_vocals -1-vs--2 split applied to the now-
        # obsolete roots_{32,128}ms_lfmethod directories, not these.)
        expected_len = WINDOW_LENGTH - 1

        valid_indices = []
        roots_list = []
        n_empty = 0
        n_bad_length = 0
        n_nonfinite = 0

        for i in range(len(raw_roots)):
            r = raw_roots[i, 0]
            if r is None:
                n_empty += 1
                continue
            r = np.asarray(r).flatten()
            if len(r) == 0:
                n_empty += 1
                continue
            if len(r) != expected_len:
                n_bad_length += 1
                continue
            if not np.all(np.isfinite(r)):
                n_nonfinite += 1
                continue
            valid_indices.append(i)
            roots_list.append(r)

        if not roots_list:
            return None

        # MIR-1K .pv/.f0 labels are MIDI note numbers (0 = unvoiced), unlike
        # PTDB's .npy which is already Hz.
        pitch_midi = np.loadtxt(label_path).reshape(-1)
        valid_indices = np.asarray(valid_indices, dtype=np.int64)

        valid_mask = valid_indices < len(pitch_midi)
        valid_indices = valid_indices[valid_mask]
        roots_list = [r for j, r in enumerate(roots_list) if valid_mask[j]]

        if not roots_list:
            return None

        roots = torch.from_numpy(np.stack(roots_list)).to(torch.complex64)
        pitch_midi = pitch_midi[valid_indices]

        T = len(roots)
        if len(roots) != len(pitch_midi):
            T = min(len(roots), len(pitch_midi))
            roots = roots[:T]
            pitch_midi = pitch_midi[:T]

        n_invalid = n_empty + n_bad_length + n_nonfinite
        if n_invalid:
            print(
                f"[ROOT WARN] {os.path.basename(root_path)}: "
                f"kept {len(roots)}/{len(raw_roots)} frames; "
                f"dropped empty={n_empty}, bad_length={n_bad_length}, "
                f"nonfinite={n_nonfinite}"
            )
        if T == 0:
            return None

        pitch_label = torch.zeros(T, self.num_class, dtype=torch.float32)
        voice_label = torch.zeros(T, dtype=torch.float32)

        for i, midi in enumerate(pitch_midi):
            if float(midi) > 0:
                freq = 440.0 * (2.0 ** ((float(midi) - 69.0) / 12.0))
                cent = 1200.0 * np.log2(freq / 10.0)
                index = int(round((cent - CONST) / 20.0))
                if 0 <= index < self.num_class:
                    pitch_label[i, index] = 1.0
                    voice_label[i] = 1.0

        return {
            'roots': roots,
            'pitch': pitch_label,
            'voice': voice_label,
            'file': root_path,
        }

    def __getitem__(self, index):
        entry = self.index[index]
        if self.use_cache:
            return torch.load(entry, map_location='cpu')
        root_path, label_path = entry
        group = next(g for g in self.groups)  # single-group experiments only when cache is off
        item = self._parse(group, root_path, label_path)
        if item is None:
            raise RuntimeError(f"{root_path} has zero usable frames")
        if self.chunked:
            T = item['roots'].shape[0]
            if T <= self.n_steps:
                return item
            begin = (index % max(1, T // self.n_steps)) * self.n_steps
            begin = min(begin, max(0, T - self.n_steps))
            end = begin + self.n_steps
            return {
                'roots': item['roots'][begin:end],
                'pitch': item['pitch'][begin:end],
                'voice': item['voice'][begin:end],
                'file': item['file'],
            }
        return item

    def __len__(self):
        return len(self.index)

    def file_at(self, index):
        return self.root_paths[index]

    def condition_at(self, index):
        return self.condition[index]


# ── PTDB SNR helpers ───────────────────────────────────────────────────────
# Support both filename conventions seen in the PTDB noisy-root generation:
#   snr00_F01_sa1.mat
#   snr00dB_F01_sa1.mat
# Clean labels remain under the original clean stem, e.g. F01_sa1.npy.
_SNR_PREFIX_PATTERNS = (
    re.compile(r"^snr(\d{2})_", re.IGNORECASE),
    re.compile(r"^snr(\d{2})dB_", re.IGNORECASE),
)


def strip_snr_prefix(stem):
    """Return (clean_stem, snr_db). Clean files return (stem, None)."""
    for pat in _SNR_PREFIX_PATTERNS:
        m = pat.match(stem)
        if m:
            return stem[m.end():], int(m.group(1))
    return stem, None


def parse_condition(file_path):
    """Return 'clean' or 'snrNN' from a PTDB root/cache source filename."""
    stem = os.path.splitext(os.path.basename(file_path))[0]
    _, snr_db = strip_snr_prefix(stem)
    return "clean" if snr_db is None else f"snr{snr_db:02d}"


class PTDB(Dataset):
    """
    PTDB-TUG dataset for the optimized RRCGD-RMVPE experiment.

    Changes relative to the old optimized loader:
      1. Uses roots_WL_<window>ms_HOP_<hop>ms, so the existing
         PTDB_data_10ms_hop/roots_WL_32ms_HOP_10ms directory is selected.
      2. clean-only vs clean+noisy is controlled by EXP_CLEAN_ONLY.
      3. Noisy roots reuse the clean .npy label by stripping the SNR prefix.
      4. .mat parsing is cached to .pt files. Training sequences are cached as
         individual fixed-length chunks so __getitem__ never re-parses .mat.
      5. Validation/test items are cached as complete files and can therefore
         be evaluated file-by-file without loading the whole dataset into RAM.

    The actual RRCGD computation remains in src/spec.py and is NOT cached here.
    The cache stores roots + pitch/voicing targets only, which removes the
    repeated MATLAB parsing / label construction CPU load without changing the
    optimized RRCGD representation.
    """

    def __init__(self, path, hop_length=None, sequence_length=None, groups=None,
                 use_cache=None):
        self.path = path
        hop_ms = HOP_MS if hop_length is None else int(hop_length)
        self.HOP_LENGTH = int(round(hop_ms / 1000 * SAMPLE_RATE))
        self.hop_ms = hop_ms
        self.seq_len = None if not sequence_length else int(sequence_length * SAMPLE_RATE)
        self.num_class = N_CLASS
        self.groups = list(groups or [])
        self.use_cache = USE_CACHE if use_cache is None else bool(use_cache)
        self.chunked = self.seq_len is not None

        if self.chunked:
            self.n_steps = self.seq_len // self.HOP_LENGTH + 1
        else:
            self.n_steps = None

        print(
            f"Indexing {len(self.groups)} group{'s' if len(self.groups) != 1 else ''} "
            f"of PTDB at {path} | window={WINDOW_LENGTH} samples "
            f"({WINDOW_LENGTH / SAMPLE_RATE * 1000:.0f} ms), hop={self.hop_ms} ms "
            f"({self.HOP_LENGTH} samples), seq_frames={self.n_steps}, "
            f"CLEAN_ONLY={CLEAN_ONLY}, cache={'on' if self.use_cache else 'off'}"
        )

        # Each entry is either a chunk-cache path (training) or a full-file
        # cache path (validation/test). root_paths keeps the original source
        # path for condition/file reporting without opening the cache.
        self.index = []
        self.root_paths = []
        self.condition = []
        skipped = 0

        for group in self.groups:
            pairs = self.files(group)
            for root_path, label_path in tqdm(
                pairs, desc=f"Indexing group {group}"
            ):
                if self.use_cache:
                    cache_entries = self._ensure_cache(group, root_path, label_path)
                    if not cache_entries:
                        skipped += 1
                        continue
                    for cp in cache_entries:
                        self.index.append(cp)
                        self.root_paths.append(root_path)
                        self.condition.append(parse_condition(root_path))
                else:
                    # Cache-off is intentionally lazy: no .mat parsing here.
                    self.index.append((root_path, label_path))
                    self.root_paths.append(root_path)
                    self.condition.append(parse_condition(root_path))

        if skipped:
            print(f"[INFO] Skipped {skipped} file(s) with no usable frames.")
        print(f"Indexed {len(self.index)} dataset item(s) for {self.groups}.")

    @staticmethod
    def availabe_groups():
        return ['train', 'valid', 'test']

    def files(self, group):
        window_ms = int(round(WINDOW_LENGTH / SAMPLE_RATE * 1000))
        roots_dir = os.path.join(
            self.path, group,
            f"roots_WL_{window_ms}ms_HOP_{self.hop_ms}ms"
        )
        if not os.path.isdir(roots_dir):
            raise FileNotFoundError(
                f"Expected PTDB roots directory not found:\n  {roots_dir}\n"
                f"Check SAMPLE_RATE={SAMPLE_RATE}, WINDOW_LENGTH={WINDOW_LENGTH}, "
                f"HOP_MS={self.hop_ms}."
            )

        root_files = sorted(glob(os.path.join(roots_dir, '*.mat')))
        if CLEAN_ONLY:
            root_files = [
                rf for rf in root_files
                if strip_snr_prefix(os.path.splitext(os.path.basename(rf))[0])[1] is None
            ]

        pairs = []
        missing = []
        for root_path in root_files:
            stem = os.path.splitext(os.path.basename(root_path))[0]
            clean_stem, _ = strip_snr_prefix(stem)
            label_path = os.path.join(self.path, group, clean_stem + '.npy')
            if not os.path.isfile(label_path):
                missing.append((stem, label_path))
            else:
                pairs.append((root_path, label_path))

        if missing:
            preview = "\n".join(f"  {s} -> {p}" for s, p in missing[:5])
            raise FileNotFoundError(
                f"{len(missing)} PTDB root file(s) have no matching clean label.\n{preview}"
            )
        return pairs

    def _cache_dir(self, group):
        window_ms = int(round(WINDOW_LENGTH / SAMPLE_RATE * 1000))
        return os.path.join(
            self.path, group,
            f"_frame_cache_WL{window_ms}ms_HOP{self.hop_ms}ms_FREF10.0_NCLASS{N_CLASS}"
        )

    def _cache_stem(self, root_path):
        return os.path.splitext(os.path.basename(root_path))[0]

    def _ensure_cache(self, group, root_path, label_path):
        """Create cache once and return a list of cache files for this item."""
        cache_dir = self._cache_dir(group)
        os.makedirs(cache_dir, exist_ok=True)
        stem = self._cache_stem(root_path)

        if self.chunked:
            # A small marker lets us know the source has already been chunked.
            manifest = os.path.join(cache_dir, stem + ".manifest")
            if os.path.isfile(manifest):
                entries = []
                with open(manifest) as f:
                    for line in f:
                        cp = line.strip()
                        if cp and os.path.isfile(cp):
                            entries.append(cp)
                if entries:
                    return entries

            item = self._parse(root_path, label_path)
            if item is None:
                return []

            T = item['roots'].shape[0]
            entries = []
            n_full = T // self.n_steps
            for i in range(n_full):
                begin = i * self.n_steps
                end = begin + self.n_steps
                cp = os.path.join(cache_dir, f"{stem}__chunk{i:04d}.pt")
                if not os.path.isfile(cp):
                    _atomic_torch_save({
                        'roots': item['roots'][begin:end],
                        'pitch': item['pitch'][begin:end],
                        'voice': item['voice'][begin:end],
                        'file': root_path,
                    }, cp)
                entries.append(cp)

            # Preserve the old optimized loader's final-remainder behavior,
            # but avoid duplicating an exactly divisible last chunk.
            if T % self.n_steps != 0:
                cp = os.path.join(cache_dir, f"{stem}__chunk{n_full:04d}.pt")
                if not os.path.isfile(cp):
                    _atomic_torch_save({
                        'roots': item['roots'][-self.n_steps:],
                        'pitch': item['pitch'][-self.n_steps:],
                        'voice': item['voice'][-self.n_steps:],
                        'file': root_path,
                    }, cp)
                entries.append(cp)

            _atomic_write_lines(manifest, entries)
            del item
            return entries

        # Full-file cache for validation/test.
        cp = os.path.join(cache_dir, stem + ".pt")
        if os.path.isfile(cp):
            return [cp]
        item = self._parse(root_path, label_path)
        if item is None:
            return []
        _atomic_torch_save(item, cp)
        del item
        return [cp]

    def _parse(self, root_path, label_path):
        """Parse one .mat + clean .npy pair into CPU tensors."""
        mat_data = scipy.io.loadmat(root_path)
        if 'roots' not in mat_data:
            raise KeyError(f"{root_path} does not contain MATLAB variable 'roots'")
        raw_roots = mat_data['roots']
        expected_len = WINDOW_LENGTH - 1

        # IMPORTANT: the MATLAB root extractor can leave individual frames
        # empty (0 roots), especially for silent/unstable frames. Do not
        # compact roots and then simply truncate the pitch labels: that would
        # shift every label after the first invalid frame. Keep the original
        # frame indices and select the corresponding labels.
        valid_indices = []
        roots_list = []
        n_empty = 0
        n_bad_length = 0
        n_nonfinite = 0

        for i in range(len(raw_roots)):
            r = raw_roots[i, 0]
            if r is None:
                n_empty += 1
                continue
            r = np.asarray(r).flatten()
            if len(r) == 0:
                n_empty += 1
                continue
            if len(r) != expected_len:
                n_bad_length += 1
                continue
            if not np.all(np.isfinite(r)):
                n_nonfinite += 1
                continue
            valid_indices.append(i)
            roots_list.append(r)

        if not roots_list:
            return None

        pitch = np.asarray(np.load(label_path), dtype=np.float64).reshape(-1)
        valid_indices = np.asarray(valid_indices, dtype=np.int64)

        # The root extractor and label generator can differ by one frame at
        # the boundaries. Keep only root indices for which a label exists.
        valid_mask = valid_indices < len(pitch)
        valid_indices = valid_indices[valid_mask]
        roots_list = [r for j, r in enumerate(roots_list) if valid_mask[j]]

        if not roots_list:
            return None

        roots = torch.from_numpy(np.stack(roots_list)).to(torch.complex64)
        pitch = pitch[valid_indices]

        T = len(roots)
        if len(roots) != len(pitch):
            T = min(len(roots), len(pitch))
            roots = roots[:T]
            pitch = pitch[:T]

        n_invalid = n_empty + n_bad_length + n_nonfinite
        if n_invalid:
            print(
                f"[ROOT WARN] {os.path.basename(root_path)}: "
                f"kept {len(roots)}/{len(raw_roots)} frames; "
                f"dropped empty={n_empty}, bad_length={n_bad_length}, "
                f"nonfinite={n_nonfinite}"
            )
        if T == 0:
            return None

        pitch_label = torch.zeros(T, self.num_class, dtype=torch.float32)
        voice_label = torch.zeros(T, dtype=torch.float32)

        for i, hz in enumerate(pitch):
            if float(hz) > 0:
                cent = 1200.0 * np.log2(float(hz) / 10.0)
                index = int(round((cent - CONST) / 20.0))
                if 0 <= index < self.num_class:
                    pitch_label[i, index] = 1.0
                    voice_label[i] = 1.0

        return {
            'roots': roots,
            'pitch': pitch_label,
            'voice': voice_label,
            'file': root_path,
        }

    def __getitem__(self, index):
        entry = self.index[index]
        if self.use_cache:
            return torch.load(entry, map_location='cpu')
        root_path, label_path = entry
        item = self._parse(root_path, label_path)
        if item is None:
            raise RuntimeError(f"{root_path} has zero usable frames")
        if self.chunked:
            # Cache-off mode still needs fixed-size sequences for DataLoader.
            T = item['roots'].shape[0]
            if T <= self.n_steps:
                return item
            begin = (index % max(1, T // self.n_steps)) * self.n_steps
            begin = min(begin, max(0, T - self.n_steps))
            end = begin + self.n_steps
            return {
                'roots': item['roots'][begin:end],
                'pitch': item['pitch'][begin:end],
                'voice': item['voice'][begin:end],
                'file': item['file'],
            }
        return item

    def __len__(self):
        return len(self.index)

    def file_at(self, index):
        return self.root_paths[index]

    def condition_at(self, index):
        return self.condition[index]

    def load(self, root_path, label_path):
        """Backward-compatible eager loader used by older scripts."""
        item = self._parse(root_path, label_path)
        if item is None:
            return []
        if not self.chunked:
            return [item]
        data = []
        T = item['roots'].shape[0]
        n_full = T // self.n_steps
        for i in range(n_full):
            b = i * self.n_steps
            e = b + self.n_steps
            data.append({k: item[k][b:e] for k in ('roots','pitch','voice') } | {'file': root_path})
        if T % self.n_steps:
            data.append({k: item[k][-self.n_steps:] for k in ('roots','pitch','voice')} | {'file': root_path})
        return data


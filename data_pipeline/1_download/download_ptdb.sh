#!/usr/bin/env bash
# Step 1: Download the PTDB-TUG speech database.
#
# PTDB-TUG (Pitch Tracking Database from Graz University of Technology)
# is distributed as a single ~3.9GB zip, no login/registration required,
# from the SPSC Graz database page:
#   https://www2.spsc.tugraz.at/databases/PTDB-TUG/
#
# Usage:
#   ./download_ptdb.sh [output_dir]   # default output_dir: ./PTDB_raw
set -euo pipefail

OUT_DIR="${1:-./PTDB_raw}"
BASE_URL="https://www2.spsc.tugraz.at/databases/PTDB-TUG"
ZIP_NAME="SPEECH_DATA_ZIPPED.zip"

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

echo "Downloading $ZIP_NAME (~3.9GB) ..."
wget -c "$BASE_URL/$ZIP_NAME"
wget -c "$BASE_URL/$ZIP_NAME.md5"
wget -c "$BASE_URL/$ZIP_NAME.sha1"

echo "Verifying checksum ..."
if command -v md5sum >/dev/null 2>&1; then
    expected=$(cut -d' ' -f1 "$ZIP_NAME.md5")
    actual=$(md5sum "$ZIP_NAME" | cut -d' ' -f1)
    if [ "$expected" != "$actual" ]; then
        echo "MD5 MISMATCH: expected $expected, got $actual" >&2
        exit 1
    fi
    echo "MD5 OK."
fi

echo "Extracting ..."
unzip -q "$ZIP_NAME" -d .

echo "Done. Raw PTDB-TUG data is under: $OUT_DIR"
echo "Expected top-level layout after extraction: 'SPEECH DATA/{FEMALE,MALE}/<speaker>/...'"
echo "Next: run ../2_split/split_ptdb.py --raw-dir \"$OUT_DIR\" --out dataset/PTDB_data_10ms_hop"

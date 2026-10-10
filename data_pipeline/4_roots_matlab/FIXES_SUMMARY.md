# Lindsey-Fox (WebLF21) — bugs found and fixed

## Bug 1 (critical — this is the one that will actually crash `lroots`)

**File:** `shared/innerparts.m` → renamed to `shared/innerpart.m`

The file was named `innerparts.m` (plural) but defined
`function [y,err] = innerpart(x, zin)` (singular). MATLAB resolves a
function call by the **filename** on the path, not by the name on the
`function` line. Every caller in the codebase — `laguerre3.m`,
`laguerre2.m`, `newton2.m`, `newton3.m` — calls `innerpart(...)`
(singular).

This function is only *supposed* to run when the compiled mex file
`mexFiles/innerpart.dll` isn't available. That `.dll` is 32-bit-Windows
only, so on Mac, Linux, or 64-bit Windows without a recompiled mex,
MATLAB falls through to the `.m` version — and with the original
filename, it can't find it at all:

```
Undefined function 'innerpart' for input arguments of type 'double'.
```

**Confirmed live in the critical path:**
`lroots.m` → `bigPolish.m` → `laguerre3.m` → `innerpart(...)`
This runs on every root, every grid-search pass. It's not an edge case.

**Fix:** renamed `innerparts.m` → `innerpart.m` (no code changes inside
the file were needed — the function body was correct all along).

## Bug 2 (latent — matters on MATLAB R2023b or newer, or at degree ≥ 2^18)

**Files:** `shared/rootsIni.m`, `lroots/lroots.m` (×2), `lroots/gridSearch.m`,
`lroots/ix2dmin.m`

All five parsed the MATLAB major version like this:
```matlab
ver = version;
ver = str2num(ver(1));   % takes only the FIRST CHARACTER of the version string
```
This worked while MATLAB's version string looked like `9.11.0...`
(`ver(1)` = `'9'`). Starting **R2023b**, MathWorks switched to a
calendar-based scheme (`23.2.0...`, `24.1.0...`), so `ver(1)` now
grabs `'2'` → parsed as major version 2, not the real value.

These checks gate `packMemory` (memory defragmentation for degree ≥
2^18 ≈ 262,144) and a loop-vs-vectorized branch in `ix2dmin.m`. At your
current scale (order ~1000+) this doesn't cause incorrect output — the
branches it guards only matter for very large-degree runs — but it's
wrong and will bite you later if you push to bigger polynomials on a
modern MATLAB, so I fixed it now.

**Fix:** replaced with a version-scheme-agnostic parse:
```matlab
ver = str2double(regexp(version,'^\d+','match','once'));
```
This grabs the full leading integer, so it's correct under both the
old and new versioning schemes.

## Bug 3 (blocks first run entirely) — hardcoded Windows log-file path

**File:** `shared/rootsIni.m`

`rootsIni` (called on the very first line of `lroots.m`) hardcoded
`folder = 'c:\poly\';` for its two log files (`roots.his`,
`DeflatedRadii.his`) and then actively tests that the folder is
writable, throwing:
```
<<Please edit RootsIni.m in folder shared and set the variable "folder"
to be a folder on your hard drive.>>
```
on any machine without that exact Windows path. The original docs
(`help/installation.htm`) say to hand-edit this every time you install
the package.

**Fix:** `folder` is now derived automatically from the package's own
location (a `logs/` folder created next to `WebLF21_fixed/`), created
with `mkdir` if it doesn't exist. No manual editing needed, and it
works the same on Windows/Mac/Linux.

Also changed `DIAGNOSTICS` from the shipped default of `1` to `0` —
with it on, every single frame in a batch run prints multi-line
progress messages, which floods the console across hundreds of frames.
The benchmark script's own progress prints are enough; flip it back to
`1` in `shared/rootsIni.m` if you want to debug one factorization in
detail.

## Critical performance note — mex compilation is not optional at your scale

`help/installation.htm` documents the package's own reference timings
(2.66 GHz single core): with the mex files compiled, degree-1000 takes
0.51s; without them (pure-.m fallback), the same degree takes 35.67s —
and the gap widens sharply with degree. Extrapolating to your actual
audio's polynomial degrees (deg 1535 at 32ms, deg 6143 at 128ms,
48 kHz, ~800 frames each) puts a whole-file `lroots` run at roughly:

|                | with mex compiled | pure-.m fallback (no mex) |
|----------------|-------------------|----------------------------|
| 32ms window    | ~9 min            | ~18 **hours** |
| 128ms window   | ~26 min           | ~245 **hours** |

(Rough power-law extrapolation beyond the package's tested range —
treat as order-of-magnitude, not exact. `benchmark_lroots_vs_roots.m`
measures your actual machine's per-frame time on a pilot frame before
committing to the full run, so you'll see real numbers, not just this
estimate.)

**Compile the mex files before running the benchmark on whole audio
files.** See `setup_lf_path.m` for the exact commands.

## Bug 4 (blocks first run entirely) — multroot toolbox missing from path

**File:** `setup_lf_path.m`

`lroots.m` unconditionally calls `extremeRoots()` (in `shared/`), which
calls `multroot()` — but `multroot()` lives in `multroot/multroot/`, a
folder `setup_lf_path.m` didn't add. Result:
```
Unrecognized function or variable 'multroot'.
Error in extremeRoots (line 80)
Error in lroots (line 135)
```
This isn't an edge case — it's on every single `lroots` call.

**Fix:** added `multroot/multroot` to `setup_lf_path.m`.

I also ran a full call-graph trace from `lroots.m` (following every
function call transitively) to check for any other missing-folder gaps
of this kind before you hit them one at a time. Confirmed `lroots.m`
only reaches `lroots/`, `shared/`, and `multroot/multroot/` — all now
on the path. (`Unwrap/`, `util/`, `broots/`, `simplified/` are for
other separate entry points — `unwrapz`, `nroots`, `broots`, the demo
script — not needed by `lroots` itself, but harmless to leave on the
path.)

## Bug 5 (crashes lroots on essentially every input) — nextpow2 called on the polynomial vector instead of its length

**Files:** `shared/iterativeDeflation.m` (line 96), `lroots/lroots.m` (line 221)

Both had `2^nextpow2(x)` where `x` (or `xx`) is the full polynomial
coefficient vector, not a scalar length. `nextpow2` applied to a
vector returns a vector (elementwise), so `2^vector` then hits
MATLAB's matrix-power operator on a non-square, non-scalar argument:
```
Error using ^
Incorrect dimensions for raising a matrix to a power. ...
```
This isn't an edge case either — `iterativeDeflation` runs on every
`lroots` call, so this reproduced on literally the simplest possible
test case, `lroots([1 -3 3 -1])`.

**Fix:** both now correctly use the *length* of the vector, matching
how every other correct `nextpow2` call site in the package does it
(e.g. `unfactor.m`: `nextpow2(degree + 1)`, `broots.m`:
`nextpow2(length(xx))`):
- `iterativeDeflation.m`: `Lsp = 2^nextpow2(x)` → `2^nextpow2(Lx)` (`Lx = length(x)` was already computed one line above)
- `lroots.m`: `nextpow2(xx)` → `nextpow2(length(xx))`

I checked every other `nextpow2(...)` call in the package (13 total)
and confirmed the rest all correctly pass a scalar (a length, a
degree, or an already-scalar variable) — these two were the only
instances of this bug.

## Bug 6 (root miscounting — traced from your own debugging session) — exact float equality on imag(z)

**Files:** `shared/hpdegree.m`, `shared/hp2fp.m`, and the final mirroring
step in `lroots/lroots.m`

Your debugging trace (frame=256, expected 255 roots, got 254) pointed
correctly at `iterativeDeflation.m`'s stopping logic. Following it
further: `iterativeDeflation` decides it's done by comparing
`hpdegree(zFar) + hpdegree(zUC)` (etc.) against the polynomial's true
`degree`. `hpdegree.m` was counting a root as "real" (contributes 1 to
the degree) only if `imag(z) == 0` **exactly** — otherwise it assumed
the root was complex with an unreturned conjugate (contributes 2).

A numerically polished root that is mathematically real essentially
never comes back with `imag()` exactly `0.0` — it carries tiny
floating-point residue (this package's own polishing routines,
`laguerre3.m`/`laguerre2.m`/`newton2.m`/`newton3.m`, and even the
third-party `multroot`'s `clipSmallParts.m`, all use a numerical
tolerance for exactly this reason — `hpdegree.m` was the one outlier
using exact equality). So a genuine real root gets counted as 2
half-degree units instead of 1, `hpdeg` reaches `degree` one root too
early, and deflation stops short. This reproduces on any polynomial
with an odd number of real roots — guaranteed for an odd-degree real
polynomial (e.g. your 256-sample / degree-255 frame).

The same exact-equality pattern existed in **`hp2fp.m`** (called
repeatedly inside `unfactor.m`, which runs every deflation iteration)
and in the final upper→full-plane mirroring step at the end of
**`lroots.m`** — fixing `hpdegree.m` alone would have just moved the
symptom: with degree-counting fixed but classification still exact,
those two would start fabricating a spurious near-duplicate conjugate
for every near-real root instead.

**Fix:** all three now classify "real" using a numerical tolerance
(`abs(imag(z)) <= 1e-8 * max(1, abs(z))`, consistent with the
tolerances the rest of the package already uses) and snap qualifying
roots to exactly real before use. I swept the rest of the codebase
(`grep` for `imag(...)==0` and `imag(...)>0`-style exact tests) and
confirmed these were the only three live sites in the `lroots` call
path — the rest were comments in the unrelated `multroot` testsuite.

Also folded in the `mexFiles/deflate.c` fix from your own compilation
session — it was missing `#include <string.h>` and `#include <math.h>`
for `memcpy`/`fabs`, which modern compilers warn on. Checked the other
four mex `.c` sources (`rdft.c`, `rdftNoTable.c`, `innerpart.c`,
`winding.c`, `ix2dmin.c`) for the same gap — they already had what they
needed.

## Not a bug, but you need this: path setup

Functions are split across `lroots/`, `broots/`, `shared/`, `Unwrap/`,
`util/`, `simplified/`. All of them need to be on your MATLAB path
before you call `lroots`. I added **`setup_lf_path.m`** — run it once
per session:
```matlab
run('setup_lf_path.m')
[rts, rtsErr, rc, outParms] = lroots(yourPolyCoeffs);
```

It also explains the mex situation: the `.dll` files in `mexFiles/`
are Windows-only. On Mac/Linux they're simply invisible to MATLAB's
path lookup (wrong extension for that platform), so it automatically
uses the pure-`.m` fallback — correct, just slower. At your degree
range (~1000+, not the 10^5–10^6 regime this package was benchmarked
at) the pure-`.m` path should be plenty fast for research use. If you
scale up later and need the mex speed on Mac/Linux, recompile from the
`.c` sources in `mexFiles/` with `mex <file>.c` (needs `mex -setup`
done first).

## What I checked and ruled out (false alarms)

I also scanned every `.m` file for the same filename/function-name
mismatch pattern that caused Bug 1. A few more turned up
(`permuteVDC.m` internally named `permuteVCD`, `expfitJWF.m` internally
named `expfit`, etc.) but in every other case the actual callers use
the **filename**, not the internal name, so MATLAB resolves them fine.
`innerpart`/`innerparts` was the only one where callers used the
internal name instead of the filename.

I did not touch the third-party `multroot/` toolbox (used for
polynomials with clustered/multiple roots) — it's not in the `lroots`
call path and wasn't part of what you asked about.

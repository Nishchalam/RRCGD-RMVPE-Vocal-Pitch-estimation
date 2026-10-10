% setup_lf_path.m
%
% Run this once per MATLAB session (or add to startup.m) before calling
% lroots(), broots(), unwrapz(), etc.
%
% WHY THIS IS NEEDED:
%   The Lindsey-Fox package spreads its functions across several
%   subfolders (lroots/, shared/, Unwrap/, util/, broots/, multroot/).
%   Functions in one folder call functions defined in the others, so
%   ALL of them must be on the MATLAB path, not just the top folder.
%
% ABOUT mexFiles/:
%   mexFiles/*.dll are compiled for 32-bit Windows only. On Mac/Linux,
%   or 64-bit Windows without recompilation, MATLAB will not find them
%   under those names and will silently fall back to the equivalent
%   pure-MATLAB .m implementations in shared/ and lroots/ (rdft.m,
%   ix2dmin.m, deflate.m, winding.m, innerpart.m, bitReverse.m).
%   This is CORRECT and requires no code change, just slower.
%
%   IMPORTANT: at degrees in the low thousands and up (typical for
%   whole-frame RRCGD polynomials), the pure-.m fallback can be VERY
%   slow -- the package's own documented benchmarks (help/installation.htm)
%   show degree-1000 taking 35.67s per call without mex vs 0.51s with mex
%   compiled (2.66 GHz reference machine). That gap only grows with degree.
%   Compiling mex is effectively required for whole-file benchmarking at
%   window lengths of tens of ms or more. Official build commands
%   (from help/installation.htm), run from inside mexFiles/:
%       mex rdft.c rdftNoTable.c
%       mex deflate.c
%       mex innerpart.c
%       mex winding.c
%   (mex must be configured with a C compiler first: run `mex -setup` once.)
%   ix2dmin.c is not in the package's official compile list -- its .m
%   fallback is apparently not the bottleneck, so it's optional.

thisDir = fileparts(mfilename('fullpath'));

addpath(fullfile(thisDir, 'lroots'));
addpath(fullfile(thisDir, 'broots'));
addpath(fullfile(thisDir, 'shared'));
addpath(fullfile(thisDir, 'multroot', 'multroot')); % required: lroots -> extremeRoots -> multroot on every call
addpath(fullfile(thisDir, 'Unwrap'));
addpath(fullfile(thisDir, 'util'));
addpath(fullfile(thisDir, 'simplified'));
% Only add mexFiles/ if you are on Windows and want the compiled .dll
% versions to take priority over the .m fallbacks for speed.
if ispc
    addpath(fullfile(thisDir, 'mexFiles'));
end

% fprintf('Lindsey-Fox package path set up from: %s\n', thisDir);
% fprintf('Run "which innerpart" to confirm the fix: it should resolve to shared/innerpart.m\n');
% fprintf('(or to mexFiles/innerpart.dll if you are on Windows and added mexFiles/ to the path).\n');

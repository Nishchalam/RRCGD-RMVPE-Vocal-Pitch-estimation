% function deg = hpdegree(z);
%
% PURPOSE:
%   Find the degree of a real polynomial given only its upper half-plane roots.
%
% REQUIRED INPUT:
%   z = The upper half-plane roots.
%
% OUTPUT:
%   deg = The degree of the real polynomial containing z.

% USES:
%   Nothing.
%
% AUTHOR:
%   James W. Fox 4/22/2001

function deg = hpdegree(z);

if(isempty(z))
   deg = 0;
   return;
end

j = find(imag(z) < 0);
if(length(j))
  error('Some root is in the lower half plane');
end

% ORIGINAL CODE used an exact test, imag(z)==0, to decide whether a root
% is real (contributes 1 to the degree) or genuinely complex (contributes
% 2, assuming its unreturned conjugate exists in the lower half plane).
% A numerically-polished root that is mathematically real will almost
% never come back with imag() exactly 0.0 -- it carries tiny floating-
% point residue (typically ~1e-13 to 1e-10 for this package's polishing
% routines; see laguerre3.m, laguerre2.m, newton2.m, newton3.m, and
% multroot's clipSmallParts.m, which all use a numerical tolerance for
% exactly this reason). With the exact test, such a root gets counted as
% complex-with-a-missing-conjugate, hpdeg reaches `degree` one root too
% early, and iterativeDeflation() stops one root short of the true
% factorization -- this reproduces reliably on odd-degree real
% polynomials (which always have an odd number of real roots).
% FIXED: use a small relative tolerance, consistent with the rest of the
% package's convention for "close enough to real".
REALTOL = 1e-8;
isReal = abs(imag(z)) <= REALTOL * max(1, abs(z));
deg = 2*length(z) - length(find(isReal));



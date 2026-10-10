% function [z, zErr] = hp2fp(z, zErr);
%
% PURPOSE:
%   Take a set of roots in the upper half-plane and supply
%   complex conjugates of the non-real roots.
%   Any roots in the lower half-plane are ignored.
%
% REQUIRED INPUT:
%   z = Roots in the upper half-plane.
%
% OPTIONAL INPUT:
%   zErr = Error estimates for z.
%
% OUTPUT:
%   z = Roots in both half-planes.
%
% COMMENTS:
%   The name means "half-plane to full plane."  This routine 
%   provides the complete set of roots of a real polynomial
%   given only the upper half-plane roots. 
%
%   For real polynomials, lroots only searches for roots in 
%   the upper half-plane. The roots in the lower half plane 
%   are the complex conjugates of the non-real roots in the
%   upper half-plane.

% USES:
%   Nothing.
%
% AUTHOR:
%   James W. Fox 5/08/2001

function [z, zErr] = hp2fp(z, zErr);

if(isempty(z))
  zErr = [];
  return;
elseif(nargout==2 & nargin==1)
  error('You must supply $2 if you want 2 output arguments.')
end

z = z(:);
if(nargin >= 2)
   zErr = zErr(:);
end

% FIXED: partition using a numerical tolerance instead of exact equality
% on imag(z). See hpdegree.m for the full explanation -- a numerically
% polished real root almost never comes back with imag() exactly 0.0,
% so the original exact test wrongly classified such roots as complex
% and fabricated a spurious near-duplicate conjugate for each one. This
% function is called from unfactor.m on every deflation iteration
% (iterativeDeflation.m), so the bug compounded across the whole run.
%
% Preserves the original contract: roots strictly in the lower half
% plane (imag < 0, and not just near-real residue) are still ignored.
REALTOL = 1e-8;
isNearReal = abs(imag(z)) <= REALTOL * max(1, abs(z));

ixImag = find(imag(z) > 0 & ~isNearReal);
complexRoots = z(ixImag);
if(nargout >= 2)
   complexRootsErr = zErr(ixImag);
end

ixReal = find(isNearReal);
realRoots = real(z(ixReal)); % snap to exactly real now that we've classified it as such
if(nargout >= 2)
   realRootsErr = zErr(ixReal);
end

Lcomplex = 2*length(complexRoots); %  2* for conjugates
Lreal = length(realRoots);
z = zeros(Lreal + Lcomplex, 1);
z(1:2:Lcomplex-1) = complexRoots;
z(2:2:Lcomplex)   = conj(complexRoots);
z(Lcomplex+1:end) = realRoots;

if(nargout >= 2)
   zErr = [complexRootsErr; complexRootsErr; realRootsErr];
end

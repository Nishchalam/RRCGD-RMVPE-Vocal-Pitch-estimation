% function [z, zErr] = nonPolish(z, zErr, conjFlag);
%
% PURPOSE:
%   Remove values from z than can not possibly be polished roots.
%   If a root has a negative imaginary part, take complex conjugate.
%
%   z is supposed to be an array of polished roots. We removed
%   trailing zeros from the coefficient array so 0 can not be
%   a root. Also, some values may have failed to polish or 
%   polished to the lower half-plane. Deal with them. 
%
%   After calling this routine, it is safe to invert the radius
%   a second time if you polished a root outside the unit circle 
%   by first flipping the coefficients and inverting the radius.
%
% REQUIRED INPUT:
%   z = Putative roots.
%
% OPTIONAL INPUT:
%   zErr - An array or error estimates for z. If [] is supplied, an array of 1e-7's
%          is used. If only a single number is supplied, it is repeated. 
%   conjFlag - Default = 1. If 1, then take the complex conjugate of any roots in
%              the lower half plane; i.e. force them into the upper half plane.
%              If not 1, do nothing.
%
% REQUIRED OUTPUT:
%   z = Believable roots. Impossible roots have been removed.
%
% OPTIONAL OUTPUT:
%   zErr = The error estimates of the believable roots.

% USES:
%   Nothing.
%
% AUTHOR:
%   James W. Fox 4/22/2001

function [z, zErr] = nonPolish(z, zErr, conjFlag);

if(nargin < 1)
   error('Supply at least 1 argument.');
elseif(nargin == 1 & nargout > 1)
   error('You must supply two input arguments if you want two output arguments.');
elseif(isempty(zErr))
   zErr = 1e-7 + zeros(size(z));
elseif(length(zErr == 1))
   zErr = abs(zErr) + zeros(size(z));
elseif(length(zErr) ~= length(z))
    error('length($3) must equal length($4).')
elseif (min(zErr) < 0)
    error('$4 must be not contain negative values.')
end

if(nargin < 3)
   conjFlag = 1; % Set default.
end

absZ = abs(z);

k = find(0<absZ & absZ<Inf & ~isnan(z));
z = z(k);
if(nargout > 1)
   zErr = zErr(k);
end

if(conjFlag == 1)
   k = find(imag(z) < 0);
   if(~isempty(k))
      z(k) = conj(z(k));
%      fprintf('          -------- FYI nonPolish found a root with negative imaginary part. Index=%g\n',k);
%      for n=1:length(k)
%         hwrite(['          -------- FYI nonPolish found a root with negative imaginary part. Index=',num2str(k(n))]);
%      end
   end
end

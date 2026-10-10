% function [z1d, z1dErr] = polishMultroot(x, z2d, ferr);
%
% PURPOSE
%   Polish any roots output by multroot that are not multiple roots.
%   (It is impossible to directly polish multiple roots accurately.)
%
% REQUIRED INPUT
%   x - any polynomial.
%   z2d - a 2d array of roots and multiplicities as output my multroot.
%   ferr - multroot forward error estimate.
%
% REQUIRED OUTPUT
%   z1d - an array of roots where a root with multiplicity k occurs k times.
%   z1dErr - error estimates for the roots

% AUTHOR
%   James W. Fox 06/05/04

function [z1d, z1dErr, relCorr] = polishMultroot(x, z2d, ferr);

if(nargin ~= 3)
   error('Supply exactly 3 arguments - the polynomial and the 2d root array from multroot.');
elseif(~isreal(x))
   error('Currently, the polynomial must have real coefficients.');
elseif(isempty(z2d))
   z1d = [];
   z1dErr = [];
   relCorr = [];
   return;
end

for k = 1:size(z2d, 1)
   if(z2d(k,2) == 1)
      [rts, zE, rC] = laguerre3(x,z2d(k,1));
      z2d(k,1) = rts;
      z2d(k,3) = zE;
      z2d(k,4) = rC;
   else
      z2d(k,3) = ferr;
      z2d(k,4) = ferr;
   end
end

[z1d, z1dErr, relCorr] = z2d_z1d(z2d);
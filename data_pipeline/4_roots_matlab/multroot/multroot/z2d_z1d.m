% function [z1d, zErr, relCorr] = z2d_1d(z2d)
%
% PURPOSE
%   Convert a 2d array of roots to a 1d array. The 2d array has unique roots in 
%   column 1 and their multiplicities in column 2.
%
% INPUT ARGUMENTS
%   z2d a 2d array of roots and multiplicities. Column 1 has the roots while column 2
%   is the multiplicities. Optionally, polishMultroot may have added two more columns
%   of root error estimates.
%
% REQUIRED OUTPUT
%   z1d a 1d array of roots where a root with multiplicity k occurs k times.
%
% OPTIONAL OUTPUT
%   zErr - Error estimates for each root. For multiple roots this is the multroot
%          forward error. For non-multiple roots this is zErr from laguerre3. These 
%          values must have been preloaded into z2d using polishMultroot.
%   relCorr - Error estimates for each root. For multiple roots this is the multroot
%          forward error. For non-multiple roots this is zErr from laguerre3. These 
%          values must have been preloaded into z2d using polishMultroot.

% AUTHOR
%   James W. Fox 04/15/04. Revised 06/05/04 to allow for 2 extra columns of error
%   estimates.

function [z1d, zErr, relCorr] = z2d_z1d(z2d)

[m,n] = size(z2d);
if(n ~= 2 & n ~= 4)
   error('%1 must have 2 or 4 columns.');
elseif(nargout > 1 & n ~= 4)
   error('$1 must have 4 columns if you want more than one output argument.');
end

rts = z2d(:,1);
multiplicities = z2d(:,2);
degree = sum(multiplicities);
z1d = zeros(degree,1);
if(n == 4)
   zErr = z1d;
   relCorr = z1d;
end

start = 1;
for k = 1:length(rts)
   mult = multiplicities(k);
   z1d(start : start+mult-1) = rts(k);
   if(n == 4)
      zErr(start : start+mult-1) = z2d(k,3);
      relCorr(start : start+mult-1) = z2d(k,4);
   end
   start = start + mult;
end

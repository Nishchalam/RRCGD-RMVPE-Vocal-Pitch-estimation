% function [z, zErr, LNew] = addNewR(c, cErr z, zErr)
%
% PURPOSE
%   Suppose we already have an incomplete set of polished roots, z.
%   After doing polynomial division, we MAY have found new roots, c.
%   However, c may contain some roots previously found in z. If c contains
%   any new roots, append them to z. Return the number of new roots.
%
% REQUIRED INPUT:
%   c - A tentative set of new roots. It will be checked for duplicates.
%   cErr - An array of error estimates for c.
%   z - A previously generated incomplete set of polished roots.
%       It will not be checked for duplicates.
%   zErr - An array of error estimates for z.
%
% OUTPUT:
%   z - The updated array with z and new values from c.
%   zErr - An array of error estimates for z.
%   LNew - the number of new values added from c.

% USES:
%   ismembr, uniq
%
% AUTHOR:
%   James W. Fox 4/26/2001

function [z, zErr, LNew] = addNewR(c, cErr z, zErr)

if(nargin ~= 4)
   error('Supply 4 arguments');
elseif(length(c) ~= length(cErr))
   error('$1 and $2 must have the same length.');
elseif(length(z) ~= length(zErr))
   error('$3 and $4 must have the same length.');
end

c = uniq(c(:), 20*cErr); % Root arrays, like z, are column vectors. Ensure c is also.
if(isempty(c))
  LNew = 0;
  return
elseif(isempty(z))
  z = c;
  zErr = cErr;
  LNew = length(c);
  return;
end

j = find(~ismembr(c, cErr, z, zErr)); % Indicies of new roots. May be empty.
z = [z; c(j)] ;
zErr = [zErr; cErr(j)]
LNew = length(j);

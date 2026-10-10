% function ixlnz = ixLastNZ(x, threshold);
%
% PURPOSE:
%   Find the index of the last non-zero value. More precisely,
%   find the last k where abs(x(k)) > threshold.
%   Return 0 if x is empty or all the entries are
%   less than or equal to the thershold.
%
% REQUIRED INPUT:
%   x - any array
%
% OPTIONAL INPUT:
%   threshold = if abs(x(k)) <= threshold, consider x(k) to be zero. 
%               Default = 0.0
%
% OUTPUT:
%   ixlnz = The index of the last non-zero value of x.

% USES:
%   nothing
%
% AUTHOR:
%   James W. Fox 4/23/2001

function ixlnz = ixLastNZ(x, threshold);

if(nargin == 1) 
  Threshold = 0.0;
else
  Threshold = abs(threshold);
end

ixlnz = find(abs(x) > Threshold);

if(isempty(ixlnz)) % x=[] or x is entirely <= threshold.
  ixlnz = 0;
else
  ixlnz = ixlnz(end);
end

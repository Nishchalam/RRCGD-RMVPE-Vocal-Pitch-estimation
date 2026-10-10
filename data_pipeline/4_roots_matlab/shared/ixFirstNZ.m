% function ixfnz = ixFirstNZ(x, threshold);
%
% Purpose:
%   Find the index of the first non-zero value. More precisely,
%   find the first k where abs(x(k)) > threshold.
%   Return 0 if x is empty or all the entries are
%   less than or equal to the thershold.
%
% REQUIRED INPUT:
%   x - Any array
%
% OPTIONAL INPUT:
%   threshold = if abs(x(k)) <= threshold, consider x(k) to be zero. 
%               Default = 0.0
%
% OUTPUT:
%   ixfnz = The index of the first non-zero value of x.

% USES:
%   Nothing.
%
% AUTHOR:
%   James W. Fox 4/23/2001

function ixfnz = ixFirstNZ(x, threshold);

if(nargin == 1) 
  Threshold = 0.0;
else
  Threshold = abs(threshold);
end

ixfnz = find(abs(x) > Threshold);

if(isempty(ixfnz)) % x=[] or x is entirely  <= threshold.
  ixfnz = 0;
else
  ixfnz = ixfnz(1);
end


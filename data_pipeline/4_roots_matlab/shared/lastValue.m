% function ans = lastValue(x);
%
% PURPOSE:
%   Find the last value of x.
%
% REQUIRED INPUT:
%   x - Any array.
%
% OUTPUT:
%   ans - The last value of x.
%
% COMMENTS:
%   This is useful for finding the last value of a complicated 
%   expression. For example, lroots has 3 lines like:
%      ixBreak1 = lastValue(find(Radii >= 1-break1*inc*radialM));

% USES:
%   Nothing.
%
% AUTHOR:
%   James W. Fox 4/22/2001

function ans = lastValue(x);

ans = x(end);

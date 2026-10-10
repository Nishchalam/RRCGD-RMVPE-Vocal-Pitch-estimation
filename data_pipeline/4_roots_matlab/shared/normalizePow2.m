% function [y, normalizationFactor] = normalizePow2(x);
%
% PURPOSE:
%   Normalize an array by dividing by a power of 2 so that if minX and maxX
%   are the minimum and maximum absolute values, then maxX is approximately
%   1/minX. This is useful for scaling polynomials so that they are less 
%   likely to underflow/overflow when you evaluate them or deflate roots.
%
% REQUIRED INPUT:
%   x - Any array.
%
% OUTPUT:
%   y - The normalized array.

% USES:
%   Nothing.
%
% AUTHOR:
%   James W. Fox 5/6/2003

% COMMENT
%   Let z=1:18; w=poly(z). Then polyval(w,z)==polyval(normalizePow2(w),z)==0. 
%   However, polyval(w/max(w),z) is never 0 when it should be. Why? 
%   NormalizePow2 only affects the exponent, not the mantissa so there is 
%   no roundoff error. By contrast w/max(w) has roundoff errors in the 
%   mantissa and the polynomial is so ill-conditioned that these tiny 
%   roundoff errors cause serious problems in the polynomial evaluation.

function [y, normalizationFactor] = normalizePow2(x);

if(isempty(x))
   y = [];
   normalizationFactor = 1;
   return;
end

absX = abs(x);
minX = min(absX);
maxX = max(absX);
if(minX == maxX) % Handles x=[Inf Inf]
   y = x;
   normalizationFactor = 1;
   return;
end

% Find the minimum non-zero absolute value.
k = find(absX);
if(length(k) == length(absX)) % No values are 0.
    minX = min(absX);
else
    minX = min(absX(k));
end

% Find the maximum non-infinite value.
k = find(absX < Inf);
if(length(k) == length(absX)) % No values are Inf
    maxX = max(absX);
else
    maxX = max(absX(k));
end

desired = exp(mean([log(minX), log(maxX)]));
normalizationFactor = 2^round(log2(1/desired)); % An integer power of 2.
y = x * normalizationFactor;

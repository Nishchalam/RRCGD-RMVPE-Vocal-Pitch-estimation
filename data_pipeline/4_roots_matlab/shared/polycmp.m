% function ans = polycmp(p1,p2);
%
% PURPOSE:
%   Find the maximum error between two polynomials after normalization.
%   Polynomials are first normalized to have max absolute value 1.
%
% REQUIRED INPUT:
%   p1 - Any polynomial.
%   p2 - Any polynomial.
%
% COMMENTS:
%   Both polarities are checked to ensure we get the smallest error.
%   It does not matter if one or both is a row or column vector.

% USES:
%   Nothing.
%
% AUTHOR:
%   James W. Fox 4/30/2001

function ans = polycmp(p1,p2);

if(length(p1) ~= length(p2))
  error('$1 and $2 must have the same length.')
end

p1 = normalizeJWF(p1(:));
p2 = normalizeJWF(p2(:));

% p1 and p2 can be the same and yet have different polarity.
% Thus we do not know if p1-p2 or p1+p2 is supposed to be small.
ans = min(max(abs(p1-p2)),max(abs(p1+p2)));
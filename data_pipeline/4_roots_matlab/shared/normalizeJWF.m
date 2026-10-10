% function ans = normalizeJWF(x, val);
%
% PURPOSE:
%   Normalize an array to have maximum absolute
%   value of 1. Do not change any value's sign.
%
% REQUIRED INPUT:
%   x - Any array.
%
% OPTIONAL INPUT:
%   val - Any positive value. After normalization the
%         maximum absolute value will be val.
%
% OUTPUT:
%   ans - The normalized array.
%
% COMMENTS:
%   This is useful for compressing two lines into one.
%   For example, in lroots there is the line:
%     xx = normalize(x(ixfirstnz : ixLastNZ(x)));
%   Without normalize(), it would take 2 lines.

% USES:
%   Nothing.
%
% AUTHOR:
%   James W. Fox 4/30/2001

function ans = normalizeJWF(x, val);

maxAbs = max(abs(x));
if(~maxAbs) % Array is all zero. Can't really scale.
  ans = x;
  return;
end
scale = 1/maxAbs;

if(nargin > 1)
  if(length(val) ~= 1)
    error('$2, val, must be a single number.');
  elseif(imag(val) ~= 0)
    error('$2, val, must be real.');
  elseif(val <= 0)
    error('$2, val, must be positive.');
  end
  scale = scale*val;
end

ans = x .* scale;

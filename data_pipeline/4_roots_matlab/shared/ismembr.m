% function ismem = ismembr(x, xErr, y, yErr);
%
% PURPOSE:
%   Find which complex values of x are elements of y. xErr and yErr
%   allow the test to be fuzzy. Return an array of 0's and 1's.
%
% REQUIRED INPUT:
%   x = Values to be tested for their elemental status.
%   xErr - An array or error estimates for x. If [] is supplied, an array of 1e-7's
%          is used. If only a single number is supplied, it is repeated. 
%   y = The container array.
%   yErr - An array or error estimates for y. If [] is supplied, an array of 1e-7's
%          is used. If only a single number is supplied, it is repeated. 
%
% OUTPUT:
%   ismem - An array of 0's and 1's. 1 = "It is an element."

% USES:
%   Nothing.
%
% AUTHOR:
%   James W. Fox 4/22/2001

function ismem = ismembr(x, xErr, y, yErr);

ismem = zeros(size(x));

if(nargin ~= 4)
  error('Supply exactly 4 input arguments.');
elseif(isempty(x))
  return;
elseif(isempty(y))
  return;
end

if(isempty(xErr))
   xErr = 1e-7 + zeros(size(x));
elseif(length(xErr) == 1)
   xErr = abs(xErr) + zeros(size(x));
else                          
  if(length(xErr) ~= length(x))
    error('length($1) must equal length($2).')
  elseif (min(xErr) < 0)
    error('$2 must be not contain negative values.')
  end
end

if(isempty(yErr))
   yErr = 1e-7 + zeros(size(y));
elseif(length(yErr) == 1)
   yErr = abs(yErr) + zeros(size(y));
else                          
  if(length(yErr) ~= length(y))
    error('length($3) must equal length($4).')
  elseif (min(yErr) < 0)
    error('$2 must be not contain negative values.')
  end
end

% Supppose x(j) is an approximation, with error estimate xErr(j),
% to the correct value, z. Suppose y(k) is a different approximation
% with error estimate yErr(k). Suppose z is on the straight line between
% x(j) and y(k), and the distance to x(j) is xErr(j) and the distance to
% y(k) is yErr(k). Then |x(j)-y(k)| = xErr(j)+yErr(k). If this is true
% for any k, then consider x(k) to be a member of y.

% I do not like this. Suppose y(k) is seriously ill-conditioned so yErr(k)
% is huge. This may make false matches with lots of other values. I think
% a better test is: Is y(k) in the circle of possibility of x(j). If so, is
% that x(j) also in the circle of possibility for that y(k). Only declare 
% a match if both are in the circle of possibility of the other. That is not
% what this code does.

for j=1:length(x)
  ismem(j) = any(abs(x(j) - y) <= (xErr(j)+yErr));
end

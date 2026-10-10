% function [a,b] = expfit(x,y)
%
% PURPOSE:
%   Compute an L2 fit of the form y(x) = a*b^x,
%
% INPUT
%   x - any array of abscissas. They do not need to be sorted.
%       They must not all be equal.
%   y - the corresponding ordinates.
%
% REQUIRED OUTPUT
%   a - the scalar multiplier
%
% OPTIONAL OUTPUT
%   b - the exponential base
%
% USES
%   linfit
%
% If more values of y(x) are positive then negative, all negative 
% values are ignored. 'a' will be positive. Otherwise, all positive 
% y values will be ignored and 'a' will be negative. There must be 
% at least 2 positive or 2 negative values.

% AUTHOR
%   James W. Fox

function [a,b] = expfit(x,y)

x = reshape(x,size(y));

if(isempty(x) | isempty(y))
   error('$1 and $2 must be non-empty.')'
elseif(length(x) ~= length(y))
   error('Must have length($1) == length($2).')
elseif(nargin ~= 2)
   error('Must have exactly 2 arguments.')
elseif(min(x) == max(x)) % x is constant. It is a vertical line.
   error('$1 is constant.');
   return;
elseif(min(y) == max(y)) % y is constant. y(x)=y(1)*1^x
   a = y(1);
   b = 1;
   return
end

j = find(y > 0);
k = find(y < 0);
if(length(j) <= 1 & length(k) <= 1)
   error('There must be at least 2 positive or 2 negative values in $2');        
end

positive = (length(j) >= length(k)); % More positive than negative values?
if(positive)
  xx = x(j);
  yy = log(y(j));
else
  xx = x(k);
  yy = log(-y(k));
end

c = linfit(xx,yy);
a = exp(c(2));
if(~positive)
   a = -a;
end

if(c(1) < log(realmax))
   b = exp(c(1));
end

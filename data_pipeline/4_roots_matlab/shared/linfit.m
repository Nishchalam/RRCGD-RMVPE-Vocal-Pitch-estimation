% function coeff = linfit(x, y, LpNorm)
%
% PURPOSE
%   Least squares linear fit y as a function of x. Return [slope, y-intercept].
%
% REQUIRED INPUT
%   x - any array of abscissas. They do not need to be sorted. They must not all be equal.
%   y - the corresponding ordinates.
%
% OPTIONAL INPUT
%   LpNorm - The Lp norm to minimize. Eventually I want to allow LpNorm to be 1
%   or Inf but currently 2 is the only value that has been implemented.
%
% REQUIRED OUTPUT
%   coeff - [slope, y-intercept]
%
% USES
%   Nothing

% AUTHOR
%   James W. Fox

function coeff = linfit(x,y,LpNorm)

NelsX = length(x);
if(NelsX ~= length(y))
   error('$1 and $2 must have the same length');
end

if(nargin < 2)
   error('Supply two arguments.'); 
else % Assume LpNorm==2
   SumX  = sum(x);
   SumY  = sum(y);
   m = (NelsX*sum(x.*y)-SumX*SumY)/(NelsX*sum(x.*x)-SumX*SumX);
   b = SumY/length(y) - m*(SumX/NelsX);  % b = avg(y)-m*avg(x)
   coeff = [m, b];
end
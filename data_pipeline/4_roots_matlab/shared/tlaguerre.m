% function [eval, dist, evalN, evalL, distN, distL] = tlaguerre2(degree, radius, Npts);
%
% Compare Newton versus laguerre2 using two tests
%   1) Then number of polynomial evaluations required to converge.
%   2) The distance polished.
% For both, smaller is better. We will run the comparison
% on a set of regularly spaced points from 0 to pi on a
% circle of user specified radius.
%
% We can compute two difference arrays:
%
%   eval = #(polynomial evaluations used by Newton) - 
%          #(polynomial evaluations used by laguerre2)
%
%   dist = (distance polished by Newton) - 
%          (distance polished by laguerre2)
%
% Both arrays will be sorted. Positive values indicate places
% where laguerre2 was better than Newton. Etc. The size of the
% difference indicates how much better one method was than the
% other.
%
% I often see values in eval near 85. This almost surely 
% means Newton went the maximum 50 iterations which is 100
% polynomial evaluations. laguerre2 went 5 iterations which
% is 15 polynomial evaluations. 85 = 100-15. So while Newton
% totally failed, laguerre2 converged in 5 iterations.
%
% Places where dist==Inf indicate Newton failed to converge in
% 50 iterations. Places where dist==-Inf indicate laguerre2 failed
% to converge in 50 iterations. I have never seen an instance
% where laguerre2 failed to converge.
%
% laguerre2() and newton3() return the number of iterations.
% laguerre2() requires 3 polynomial evaluations per iteration.
% newton3() requires 2 polynomial evaluations per iteration.
% Thus, for each point,
% 	 eval = 2*#(Newton iterations) - 3*#(laguerre2 iterations)
% The number of laguerre2 iterations must be significantly smaller
% than the number of Newton iterations if this is to be positive,
% which indicates that laguerre2 was better.
%
% You can either supply a specific polynomial or let it
% generate a random coefficient polynomial. In the latter case,
% each time this is run a new polynomial is generated, even though
% the input parameters are the same. Thus it is useful to to hit the
% up-arrow key and rerun a it few more times with the same parameters.
%
% REQUIRED INPUT
%   degree - The polynomial degree or the polynomial itself.
%              If length(degree)==1, a random coefficient polynomial
%            with this degree will be created, uniformly distributed
%            between -1 and 1,
%              If length(degree)>1, this array will be used
%            as the polynomial coefficients.
%   radius - The radius of the circle to evaluate points on
%   Npts - The number of regularly spaced points on each circle
%
% REQUIRED OUTPUT
%   eval - #(polynomial evaluations used by Newton) - 
%          #(polynomial evaluations used by laguerre2)
%
% OPTIONAL OUTPUT
%   dist  - (distance polished by Newton) - (distance polished by laguerre2)
%   evalN - The number of polynomial evaluations used by Newton.
%   evalL - The number of polynomial evaluations used by laguerre2.
%   distN - The distance polished by Newton.
%   distL - The distance polished by laguerre2.
%
%   distN, distL, evalN, and evalL are not sorted. However, 
%   eval=sort(evalN-evalL) and dist=sort(distN-distL).
%
% SAMPLE INVOCATION
%   [eval, dist]=tlaguerre2(1000,.99,24) % 24 points on 1000 degree, radius=.99
%
% Using 24 means all the output is displayed on one screen. You may want to use
% additional output args to analyze the first two output args in greater detail.
%
% WRITTEN BY
%   James W. Fox 09-05-01

function [eval, dist, evalN, evalL, distN, distL] = tlaguerre2(degree, radius, Npts)

if(nargin ~= 3)
   error('Supply exactly 3 arguments.');
elseif(length(degree)==1 & degree <= 0)
   error('The polynomial degree must be > 0.');
elseif(length(radius) ~= 1)
   error('$2 must have 1 value.');
elseif(radius <= 0)
   error('$2 must be > 0.');
elseif(length(Npts) ~= 1)
   error('$3 must have 1 value.');
elseif(Npts < 1)
   error('$3 must be > 0.');
end

if(length(degree) == 1)
   x = rand(1,degree+1);
   x = 2*(x - mean(x)); % Uniformly distributed between -1 and 1.
else
   x = degree;
end
x = normalize(x);

Npts = fix(Npts);
theta = ((0:Npts-1)*pi/(Npts-1))'; % 0 to pi radians in Npts steps.
start = radius*(cos(theta) + i*sin(theta)); % start positions for polishing.

[newtZ, newtCorr, newtN] = newton3(x,start);
[lagZ, lagCorr, lagN] = laguerre2(x,start);

% Compare the number of polynomial evaluations used by the two methods.
evalN = 2*newtN'; % newton3() uses 2 polynomial evaluations per iteration.
evalL = 3*lagN';  % laguerre2() uses 3 polynomial evaluations per iteration.
eval = sort(evalN - evalL);

% Compare how far the roots polished by the two methods.
distN = abs(start-newtZ)';
distL = abs(start-lagZ)';
dist = sort(distN - distL);

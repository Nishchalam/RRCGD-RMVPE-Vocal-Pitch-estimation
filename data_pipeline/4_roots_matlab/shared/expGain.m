% function [x, exponGain, multiplier] = expGain(x)
%
% PURPOSE
%   After normalization, if the maximum of the first half is greater than 10 times
%   the maximum of the second half (or vice versa) find and apply an exponential
%   gain control. I.e., fit an exponential function to abs(x) and divide by the 
%   exponential to make the time series flatter.
%
%   Note, if a exponential gain is applied to a polynomial and the new polynomial
%   is factored to get a set of roots, z, then the roots of the original are
%   exponGain*z. Thus, to factor [1 0 0 0 -1e-100], apply expGain() gain to get
%   [1 0 0 0 -1] and exponGain=1e-25. Thus the roots of the original polynomial are
%   1e-25*(the fourth roots of unity). Well it is not exactly that; there is a 
%   roundoff error of 5 in the 16th digit.
%
% INPUT
%   x - any polynomial
%
% REQUIRED OUTPUT
%   x - the polynomial after gain control
%
% OPTIONAL OUTPUT
%   exponGain - the output polynomial was divided by exponGain.^(0:(length(x)-1))
%   multiplier - the fit is of the form multiplier*exponGain.^(0:length(x)-1)

% USES
%   expfitJWF, normalizepow2

% AUTHOR
%   James W. Fox 5/30/03

function [x, exponGain, multiplier] = expGain(x)
if(length(x) < 2)
   exponGain = 1;
   multiplier = 1;
   return;
end
x = normalizePow2(x);
Lx = length(x);
exponGain = 1; % Default in case we do not need any gain control
absX = abs(x);
max1 = max(absX(1:floor(Lx/2)));  % Maximum of the first half.
max2 = max(absX(1+floor(Lx/2):end)); % Maximum of the second half.
if(max1 > 10*max2 | max2 > 10*max1)
   t = 0 : (Lx-1);
   t = reshape(t, size(absX));
   [multiplier, exponGain] = expfitJWF(t, absX);
   if(exponGain^t(end) ~= 0)
      x = x ./ (multiplier*(exponGain.^t));
   else
      k = find(x ~= 0);
      sgn = sign(x(k));
      x(k) = 1/multiplier*real(exp(log(sgn.*x(k)) - t(k)*log(exponGain)));
   end
end

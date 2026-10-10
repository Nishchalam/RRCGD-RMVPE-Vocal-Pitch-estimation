% function [f, errF] = hornerWithError(p, z);
%
% PURPOSE
%   Use Horner's method to evaluate a polynomial at an array of values.
%   Optionally return an estimate of the error in the evaluation.
%
%   This will not be used if there is a compiled mex file for innerpart().
%
% REQUIRED INPUT ARGUMENTS
%   p - The polynomial coefficients. p(1)=highest order term.
%   z - The places where the polynomial is to be evaluated.
%
% REQUIRED OUTPUT ARGUMENTS
%   f - The polynomial evaluations.
%
% OPTIONAL OUTPUT ARGUMENTS
%   errF - Estimated error for each evaluation. 
%
% ALGORITHM
%   The calculation of the error is based on Higham, 
%   "Accuracy and Stability of Numerical Algorithms," edition 1, p. 106.
%   See Higham p. 3 for the definition of u and p. 69 for the definition of gamma2.

% AUTHOR
%   James W. Fox 6/25/2003

% By the way, it is possible to evaluate the derivative in the same loop that evaluates
% the polynomial. When tried on a random coefficient polynomial, it differed from
% polyval(polyder(p),z) in the 16'th digit. However, when tried on the Wilkinson polynomial,
% it differed in the 3'd digit! I do not know for a fact, but I bet using polyder is
% more accurate.

function [f, errF] = hornerWithError(p, z);
if(nargin ~= 2)
  error('\t<Must supply at exactly 2 arguments.>')
end

if(nargout <= 1)
   f = polyval(p,z);
   %   for j = 1:length(z) % The usual Horner evaluation
   %      y = p(1);
   %      x = z(j);
   %      for k = 2:length(p)
   %        y = x*y + p(k);
   %      end
   %      f(j) = y;
   %   end
else
   f = zeros(size(z));
   errF = f;
   u = 2^(-53);
   gamma2 = 2*u/(1-2*u);
   for j = 1:length(z)
      x = z(j);
      absX = abs(x);
      y = p(1);
      mu = abs(y)/2;
      for k = 2:length(p)
         y = x*y + p(k);
         mu = absX*mu + abs(y);
      end
      f(j) = y;
      errF(j) = sqrt(2)*gamma2*(2*mu - abs(y));
   end
end

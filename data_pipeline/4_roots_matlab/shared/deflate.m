% function yout = deflate(yin,z)
%
% PURPOSE
%   Deflate a single root of a real polynomial. If the root is
%   complex, also deflate its complex conjugate. If the root is
%   outside the unit circle, use backward deflation.
%   The output will be a real polynomial.
%
% REQUIRED INPUT
%   yin - The polynomial.
%   z   - The root to be deflated.
%
% REQUIRED OUTPUT
%   yout - The deflated polynomial.
%
% AUTHOR
%   James W. Fox

function yout = deflate(yin,z)

if(~isreal(yin))
  error('$1 must be real.');
elseif(length(z) ~= 1)
  error('$2 must have exactly 1 value.');
end

s = size(yin); % Used at the very end.
yin = yin(:).'; % Must be same structure as denom which I make a row vector.

if(isreal(z))
  denom = [1, -z];
else
  denom = [1, -2*real(z), real(z)*real(z) + imag(z)*imag(z)]; 
end

% Originally I always used deconv, not deconvFFT. It factored a 20,000 degree poly.
% However, it took 30 seconds to do 1 deconvolution. deconvFFT is 11 times faster
% to do one root of a 20,000 degree polynomial. On the 300 MHz computer, the 
% 20,000 degree took 4.2 hours using deconvFFT while deconv took 26 hours.
% However, deconv is actually faster for short yin. 

degree = length(yin) - 1;
if(abs(z) <= 1)
   if(degree < 1000)
      yout = deconv(yin, denom);
   else
      yout = deconvFFT(yin, denom);
   end
else
   if(degree < 1000)
      yout = flip(deconv(flip(yin), flip(denom)));
   else
      yout = flip(deconvFFT(flip(yin), flip(denom)));
   end
end

% deconvR is my new Mex file alternative to deconvFFT and deconv.
% It does polynomial division in the time domain.
% While deflating a real root from a 20,000 degree polynomial, deconvR
% got exactly the same answer as deconv. However, the times were:
%   deconvR    0
%   deconvFFT  1.04
%   deconv    25.92 
% deconvR differed from deconvFFT by 1.5E-13

%if(abs(z) <= 1)
%   yout = deconvR(yin, denom);
%else
%   yout = flip(deconvR(flip(yin), flip(denom)));
%end

if(s(2)==1) % yin was initially a column vector
   yout = yout.';
end


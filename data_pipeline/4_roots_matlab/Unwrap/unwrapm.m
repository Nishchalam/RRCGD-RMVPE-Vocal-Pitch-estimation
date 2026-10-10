% function [ph,pad] = unwrapm(x,pad)
%
% PURPOSE:
%   Estimate the unwrapped phase of x by repeatedly padding with zeros
%   to double the length and using Matlab's unwrap() until two
%   consecutive answers are equal at the central frequency.
%
%   This is not as accurate as unwrapz; it can be off by multiples of 2*pi.
%   However, it is easier and faster to use and even if it is wrong, it is
%   likely to be wrong by small multiples of 2*pi.
%
%   Unwrapph uses a somewhat similar strategy, but it only does one long FFT.
%   For random coefficient polynomials, if you use the defaults, unwrapm gets 
%   the correct answer 85% of the time. Unwrapph gets the correct answer 99% of 
%   the time, but it is much slower.
%
%   If there is a root on the unit circle, only unwrapz gets the correct answer.
%
% REQUIRED INPUT:
%   x - any real function
%
% OPTIONAL INPUT:
%   pad - Default=2. Round this value up to a power of 2 and zero pad to this length.
%
% OUTPUT:
%   ph - The approximate unwrapped phase
%   pad - The final pad factor that was used.
%
% USES:
%   angl, adjustUnwrapping

% AUTHOR:
%   James W. Fox. 2/15/02

% This code is freely released for research exchange only. The author 
% is not responsible for any damage from using this code.

function [ph,pad] = unwrapm(x,pad)

DIAGNOSTICS = 1; % 0=none, 1=minimal, 2=full diagnostic messages

if(~isreal(x))
   error('$1 must be real');
end

sizeX = size(x); % Used at the very end.
x = x(:)';
NminphaseOld = realmax; % Estimated # of minphase roots.
NminphaseNew = -realmax; % Quit when NminphaseOld==NminphaseNew
if(nargin == 1)
   pad = 2; % Pad with zeros to multiply the length by this factor.
else
   pad = 2^round(log2(pad));
end

while(NminphaseOld ~= NminphaseNew)
	xPadded = [x, zeros(1,length(x)*(pad-1))]; % length(xPadded)=pad*length(x)
   ph = unwrap(angl(fft(xPadded)));
   NminphaseOld = NminphaseNew;
   NminphaseNew = round(ph(length(ph)/2 + 1)/pi); % length(ph) is even - a power of 2.
   if(DIAGNOSTICS)
      fprintf(1,'Estimated # of minphase roots = %g\n',-NminphaseNew)
   end
   pad = 2*pad;
end

pad = pad/2;
ph = ph(1:pad:length(ph)); % Get a phase function that is the same length a x.

% Adjust ph so that it differs from the wrapped phase by exactly a multiple of pi.
s = fft(x);
ph = adjustUnwrapping(ph, abs(s), angl(s));
ph = reshape(ph, sizeX);

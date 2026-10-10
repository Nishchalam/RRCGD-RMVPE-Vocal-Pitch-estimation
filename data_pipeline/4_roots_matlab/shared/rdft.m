% function rdft(x, direction, nels)
%
% The mex file is a split radix real FFT without table lookup
% based on Ooura's Mathematical Software Packages
%   http://momonga.t.u-tokyo.ac.jp/~ooura/
%
% If you do not use the mex file, the code below uses Matlab's
% fft and ifft to simulate the same thing.
%
% x - The time series.
% direction =  1 - Forward real FFT
%           = -1 - Inverse real FFT
% nels - Round this up to a power of 2. Then pad x with zeros
%          to this length. Default=length(x).
%
% If direction==-1, length(x) must be 1 + a power of 2.
% However, I do not check that the first and last samples are
% real. I just ignore any imaginary parts.
% 
% If direction==-1, I do not allow you to use nels. Nels is only
% for padding in forward transforms.
%
% Warning, Ooura uses exp(+2*pi*i*omega) while Matlab uses
% exp(-2*pi*i*omega). Thus their forward transforms differ by
% complex conjugation. I could fix this in the C code for rdft,
% but that slows things down and it really does not matter.
%
% rdft is 12% faster than Matlab's fft for length 128;
% 40% faster for length 1000; 9.6 times faster for length 2^19.
%
% rdft is a more accurate transform pair than Malab's fft and ifft.
%
% AUTHOR
%  James W. Fox

%EDU» x=rand(1,2^10);
%EDU» max(abs(x-ifft(fft(x))))
%	ans =  4.7766e-015
%EDU» max(abs(x-rdft(rdft(x,1),-1)))
%	ans =  4.9960e-016

%EDU» x=rand(1,2^13);
%EDU» max(abs(x-ifft(fft(x))))
%	ans =  9.5720e-015
%EDU» max(abs(x-rdft(rdft(x,1),-1)))
%	ans =  8.0491e-016

%EDU» x=rand(1,2^16);
%EDU» max(abs(x-ifft(fft(x))))
%	ans =  1.6132e-013                Yuk!
%EDU» max(abs(x-rdft(rdft(x,1),-1)))
%	ans =  1.1102e-015 

%EDU» x=rand(1,2^20);
%EDU» max(abs(x-ifft(fft(x))))
%	ans =  5.3684e-012                Yuk!
%EDU» max(abs(x-rdft(rdft(x,1),-1)))
%	ans =  1.6896e-015

function z = rdft(x, direction, nels)

if(length(x) < 2)
   z = x;
   return;
end

s = size(x); % Used at the very last line.
x = x(:);

if(nargin < 2)
   error('Supply at least 2 arguments.');
elseif(direction ~= -1 & direction ~= 1)
   error('$2 must equal 1 or -1');
elseif(nargin==3 & nels < length(x))
   error('$3 must be >= length($1)');
elseif(direction == -1)
   LxM1 = length(x) - 1;
   if(LxM1 ~= 2^nextpow2(LxM1))
      error('length($1) must equal 1+power of 2 for inverse transform.');
   end
end

if(direction == 1)
   if(nargin < 3)
      nels = length(x);
   end
   NP2 = 2^nextpow2(nels);
   z = fft(x, NP2);
   z = conj(z(1:(NP2/2+1)));
else % direction == -1
   j = 2 : LxM1;
   z = real(ifft([x; flip(conj(x(j)))]));
end

if(s(1) == 1) % If the input was a row vector, make the output a row vector.
   z = z.';
end

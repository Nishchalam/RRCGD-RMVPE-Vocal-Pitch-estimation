% function [nEstimated, pval, lastLength] = nroots(x, maxR, startLength, stopDiff)
% 
% PURPOSE
% -------
% Given a polynomial, x, use the winding number function, winding(), to find
% the approximate number of roots whose magnitudes are <= maxR. maxR is allowed
% to be an array. 
%
% The following discussion is only for random coefficient polynomials.
%
% nroots() is extremely fast if maxR is significantly different from 1.0 but it
% is slow near 1 because most roots of random coefficient polynomials are near the
% unit circle. Near 1.0 I suggest either nroots(x, 1, 8*length(x)) or 
% nroots(x, 1, 128*length(x),Inf). I prefer the former although I have seen several
% examples where you needed 128*length(x) if you want the exactly correct answer.
%
% Let Lx=length(x). I have checked degree=1,000, 10,000, and 32,000 and nroots(x,1,128*Lx,Inf)
% almost ALWAYS gets the correct answer. However using startLength=8 will usually find the
% correct answer at 32*Lx or 64*Lx and be faster. You need to start at 8*Lx because 16*Lx
% is sometimes the correct answer but the default for stopDiff requires two consecutive equal values.
%
% If maxR<1, the answer will usually be correct or slightly too small; but not always.
% If maxR>1, the answer will usually be correct or slightly too large; but not always.
%
% REQUIRED INPUT ARGUMENTS
% ------------------------
% x - The polynomial
% maxR - An array of radii to check for the number of roots < these radii.
%
% OPTIONAL INPUT ARGUMENTS
% ------------------------
% startLength - Default=(length(x) rounded up to a power of 2). 
%     First FFT length. It will be rounded up to a power of 2 at least as big as length(x).
%     The program derives several estimates of the number of roots by successively using
%     FFT's of length: startLength, 2*startLength, 4*startLength, ... Argument stopDiff 
%     determines when it stops.
%     For some unknown reason, setting startLength to an extremely large value, i.e.,
%     2^19, can cause answers that are way too large or way too small. I suspect that
%     rdft is not sufficiently accurate for this length.
% stopDiff - Default=0. Stop when consecutive estimates differ by at most stopDiff.
%     stopDiff = 0 means stop when consecutive estimates agree exactly. This setting 
%                almost always gets the correct answer. However, it can be wrong if
%                there are roots with magnitude very very slightly smaller than MaxR.
%     stopDiff = 2 usually gets the correct answer. If it is wrong, it is usually
%                only wrong by 2. However once I have seen it be wrong by 6.
%                Eg., I have seen 3 instances where successive estimates of
%                the number of roots were like 554,556,558,560. stopDiff=2, caused
%                it to stop at 556 while stopDiff=0 got the correct answer.
%     stopDiff = Inf means stop after the first FFT. Do not do a second FFT.
%
% OUTPUT ARGUMENTS
% ----------------
% n - The number of roots with magnitude < maxR.
%     n = Inf if every number is a root; i.e. all the coefficients are 0.0.
% pval - The last set of polynomial evaluations.
% lastLength - the size of the final FFT. It will be 0 if every coefficient was 0
%     or if it was a constant polynomial. On my system, Matlab will not allow FFT's
%     of arrays of length 2^25 or larger - "Out of memory."
%
% SAMPLE USAGE
% ------------
% nroots(x,.99)       - Estimate the number of roots <= .99
% nroots(x,1,16000)   - Estimate the number of roots <= 1. Let the first FFT length be 16384.
% nroots(x,1,16000,2) - Estimate the number of roots <= 1. Let the first FFT length be 16384.
%                       Stop if two successive estimates differ by 2.
% nroots(x,1,64000,Inf) - Estimate the number of roots <= 1. Use an FFT length of 65536.
%                       Do not use a second FFT.

% USES:
%   ixFirstNZ, ixLastNZ, normalizePow2, rdft, winding
%
% AUTHOR:
%   James W. Fox 4/24/2002
%
% COMMENT ON ACCURACY
%   If r<1, then nroots(x,r) is usually correct or slightly too small. But not always.
%   Here is an odd counter example. Let coeff and rts be the data in degree200.zip. 
%   Sort the roots by magnitude. There are consecutive upper half-plane roots with
%   magnitudes 0.9877 & 0.9886. length(find(abs(rts)<.988437)==64.
%   Unfortunately nroots sometimes gets the wrong answer, even with a large zero pad
%   and sometimes it is too large. I do not know the cause.
%
% >> nroots(coeff,.9884375000000000,2^13) % Correct
%	ans =    64
% >> nroots(coeff,.9884375000000001,2^13) % It is too LARGE.
%	ans =    66
% >> nroots(coeff,.9884375000000002,2^13) % It is too small.
%	ans =    62
% >> nroots(coeff,.9884375000000002,2^14) % Redo with larger pad. It is even samller!
% 	ans =    60
% >> nroots(coeff,.9884375000000002,2^17) % Redo with larger pad. It is slow but correct.
%	ans =    64
% >> nroots(coeff,.9884375000000003,2^13) % Correct again at larger magnitude.
%	ans =    64
% >> nroots(coeff,.9884375000000004,2^13) % Correct again at larger magnitude.
%	ans =    64

function [nEstimated, pval, lastLength] = nroots(x, maxR, startLength, stopDiff)
global RAM MAX_ARRAY_LENGTH_EXPONENT

DIAGNOSTICS = 0; % 0=none, 1=minimal, 2=full diagnostic messages
rootsIni;

maxFFTLength = 2^MAX_ARRAY_LENGTH_EXPONENT;
LmaxR = length(maxR); % The number of search radii.
s = size(x);
if(nargin < 2 | nargin > 4)
   error('Supply 2, 3, or 4 arguments');
elseif(isempty(x))
   nEstimated = zeros(1,LmaxR);
   lastLength = zeros(1,LmaxR);
   return;
elseif(s(1) > 1 & s(2) > 1)
   error('$1 must be 1xn or nx1.');
elseif(~isreal(x))
   error('$1 must be real.')
elseif(min(maxR) < 0)
   error('$2 must be non-negative.');
elseif(max(abs(x)) == 0)
   nEstimated = Inf + zeros(1,LmaxR);
   lastLength = zeros(1,LmaxR);
   k = find(maxR == 0); 
   if(~isempty(k))
      nEstimated(k) = 1; % There is one root, namely 0.0.
   end
   return;
elseif(length(x) == 1) % x==0 was handled by the previous case.
   nEstimated = zeros(1,LmaxR);
   lastLength = zeros(1,LmaxR);
   return;
end   

ixlnz = ixLastNZ(x);
zeroMult = length(x) - ixlnz; % The multiplicity of 0.
x = x(ixFirstNZ(x) : ixlnz); % Remove leading & trailing zeros.
degree = length(x) - 1;
x = normalizePow2(x(:)'); % Make it a normalized row vector.

LxPow2 = 2^(1+nextpow2(length(x))); % 2*(length(x) rounded up to a power of 2).
if(LxPow2 < 16) % I have no idea why, but if x=[1 -1e-30 1e-60] and you let LxPow2
   LxPow2 = 16; % be 4, then nroots(x,.1) is 0 instead of 2 as it should be.
end

if(nargin < 3)
   startLength = LxPow2; % Default
else
   if(length(startLength) ~= 1)
     error('$3 must be a scalar.');
   end
   startLength = 2^nextpow2(max(startLength,length(x))); % Make it a power of 2.
end

if(nargin < 4)
   stopDiff = 0; % Default
else
   stopDiff = round(stopDiff);
   if(length(stopDiff) ~= 1)
     error('$4 must be a scalar.');
   elseif(stopDiff < 0)
      error('$4 must be non-negative.');
   end
end

% If maxR<1 I do the straightforward computation. If maxR>1 I first find the
% number of roots of the flipped polynomial with magnitude < 1/maxR. This is the
% same as the number of roots of the original polynomial with magnitude>maxR.
% Then I subtract this value from the degree. This has two advantages.
% 1) It is guaranteed not to overflow. If I use the straightforward computation
% when maxR>1, then I would compute maxR^degree and this could overflow.
% 2) For whatever reason it is usually 4-5 times faster; i.e. smaller values of 
% the FFT length are sufficient so fewer iterations are required.

% If you evaluate the polynomial at closely spaced values on the circle with radius 
% maxR, the winding number of this sequence of complex values should equal
% the number of roots inside the circle of radius maxR. Unfortunately, Matlab
% chose to decree that the first coefficient was the highest order coefficient.
% Because of this, the winding number is actually the number of roots outside the
% circle of radius maxR. To correct for this, simply flip() the coefficients.

lastLength = zeros(1,length(maxR)); % Allocate the 2 output arrays.
nEstimated = zeros(1,length(maxR));

for k = 1:LmaxR
  if(DIAGNOSTICS)
     disp(' ');
  end
  r = maxR(k);
  nRtsOld = -Inf; % nRtsOld = the previous estimate. Initialize to an impossible value.
  FFTLength = min(startLength, maxFFTLength);
  done = 0;
  if(r == 0)
     done = 1;
     nRts = length(x) - ixLastNZ(x); % The number of trailing zeros.
  elseif( r <= 1)
      y = (r).^(0:degree) .* flip(x); % We do flip because (1/r).^(0:degree) might overflow.
  else   
      y = (1/r).^(0:degree) .* x;
  end
  
  while(~done)
     if( r <= 1)
      pval = rdft(y,1,FFTLength);
      nRts = abs(winding(pval,0));
%     nRts = abs(winding(fft(y,FFTLength),1)); % If you insist on fft() instead of rdft().
     else
      pval = rdft(y,1,FFTLength);
      nRts = degree - abs(winding(pval,0));
%     nRts = degree - abs(winding(fft(y,FFTLength),1)); % If you insist on fft() not rdft().
    end
    
    if(DIAGNOSTICS)
       fprintf('FFT length=%6d\tEstimated # of roots=%d\n',FFTLength,nRts+zeroMult);
    end
    
    % Warning for 1,000-degree polynomials, I have seen nroots(x,1,length(x),Inf)
    % return 0, which is absurd. Do not trust nRts==0 if r is very close to 1.0.
    if(nRtsOld>-Inf & ((nRts<nRtsOld & r<1)|(nRts>nRtsOld & r>1))) % Happens if you pad too much.
      done = 1;
      if(r<1) % Choose the smaller of the two possible answers.
         nEstimated(k) = nRts + zeroMult;
         lastLength(k) = FFTLength;
      else
         nEstimated(k) = nRtsOld + zeroMult;
         lastLength(k) = FFTLength/2;
      end         
    elseif(abs(nRts-nRtsOld)<=stopDiff | (nRts==0 & abs(r-1)>5/degree) | FFTLength==maxFFTLength)
      done = 1;
      nEstimated(k) = nRts + zeroMult;
      lastLength(k) = FFTLength;
    else % Get ready to try again.
      FFTLength = 2*FFTLength;
      nRtsOld = nRts;
    end
  end
end

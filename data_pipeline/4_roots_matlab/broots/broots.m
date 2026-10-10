% function [raw, polished, err, rc, outparms] = broots(x, ID_NUM);
%
% PURPOSE
% -------
% Polynomial factorization by the Burrus argument randomization method.
%
% REQUIRED INPUT ARGUMENTS
% ------------------------
% x - the polynomial
%
% OPTIONAL INPUT ARGUMENTS
% ------------------------
% ID_NUM - Default=99. An identifying integer used by
%          iterativeDeflation() to print success/failure messages.
%
% REQUIRED OUTPUT ARGUMENTS
% -------------------------
% raw - roots obtained by deflating one root at a time but not
%       polished against the original polynomial.
%
% OPTIONAL OUTPUT ARGUMENTS
% -------------------------
% polished - raw roots after polishing against the original polynomial,
%            and then removing false duplicates, and doing a final
%            massive deflation to find missing stragglers. If you do
%            not ask for polished, no polishing will be done and the
%            job will be much quicker, but not confirmed correct.
% err - Error estimates for the polished roots.
% rc - the return code, see iterativeDeflation() for the values
% outparms - see iterativeDeflation() for the values.
%
% USES
% ----
% bitReverse, expGain, extremeRoots, hpdegree, hp2fp, iterativeDeflation,
% ixFirstNZ, ixLastNZ, laguerre2, laguerre3, normalizePow2, rootsIni, uhp, uniq

% AUTHOR
% ------
% James W. Fox 10/9/2002
%
% HISTORY
% -------
% C. Sidney Burrus provided the idea of randomizing the starting angles 
% when searching for new roots. He also provided a working Matlab program
% which consistently worked for 1,000 degree polynomials.
%
% Gary A. Sitton rewrote this program to increase the speed. It sometimes worked
% for degree 10,000. He also provided a crucial observation about deflation.
%
% James W. Fox wrote the program below. The new ideas are:
%   1) Use Laguerre instead of Newton for polishing.
%   2) Use the van der Corput, bit-reversal, method to randomize the
%      starting angles when searching for new roots.
%   3) Use recently deflated roots to determine starting magnitudes when
%      searching for the next root.
%      This version worked 9 times out of 12 for degree 100,000 but never
%      worked for degree 150,000.

% This program, has successfully factored:
%  20 polynomials of degree 10,000, 
%   3 polynomials of degree 20,000, 
%   3 polynomials of degree 50,000, 
%  14 polynomials of degree 70,000, 
%   9 polynomials of degree 100,000.
%
% It failed on:
%   1 polynomial of degree 30,000, 
%   1 polynomial of degree 70,000, 
%   3 polynomials of degree 100,000 
%   3 polynomials of degree 150,000.
%
% The one failure at 30,000 is most disappointing. It missed 2 pairs of 
% complex conjugates that were extremely close to -1. I.e., after deflating
% the roots that were found, the resulting polynomial was approximately (z+1)^4.
% The roots of this polynomial could not be found with sufficient accuracy that
% they could be polished against the original polynomial and not move far.

function [raw, polished, err, rc, outparms] = broots(x, ID_NUM);

DIAGNOSTICS = 0; % 0=none, 1=minimal, 2=full diagnostic messages

% Initialize globals
rootsIni;

% Set ID_NUM
if(nargin == 0 | nargin > 2)
   error('Supply 1 or 2 arguments.');
elseif(nargin == 1)
   ID_NUM = 0;
else
   if(length(ID_NUM) ~= 1)
      error('length($2) must be 1.');
   end
end

s = size(x);
if(s(1) ~=1 & s(2) ~= 1)
  error('$1 must be 1xn or nx1.')      
elseif(isempty(x))
  raw = [];
  polished = [];
  rc = 0;
  outParms = [0,Inf];
  return;
elseif(length(x)==1 & x~=0)
  raw = [];
  polished = [];
  rc = 0;
  outParms = [0,Inf];
  return;
elseif(~max(abs(x)))
  disp('Polynomial is entirely zero. Every number is a root.');
  raw = [];
  polished = [];
  rc = 0;
  outParms = [0,Inf];
  return;
end

if(ID_NUM)
%  fprintf('\n<<<<< broots beginning polynomial # %d >>>>>\n', ID_NUM) ;
end

%---- Remove leading and trailing zeros and normalize ----
ixLNZ = ixLastNZ(x); % Index of the last non-zero coefficient.
NZeros = length(x) - ixLNZ; % Multiplicity of 0.0 as a root. Used at the very end.

x = normalizePow2(x(ixFirstNZ(x) : ixLNZ)); % Strip leading and trailing zeros.
xx = x(:).'; % Keep x untouched, so you can polish against it near the bottom.
degree = length(xx)-1;
raw = zeros(degree,1); % Space for the answer.

% I have seen exponential gain sabotage the search for extreme roots which are sensitive.
% Therefore do not do exponential gain before extremeRoots as I originally did.

exponGain = 1; % Use this if you comment out expGain()
nFound = 0; % Use this if you comment out extremeRoots()

[extremeRts, extremeRtsErr, xx] = extremeRoots(xx);
[extremeRts, extremeRtsErr] = hp2fp(extremeRts, extremeRtsErr);
nFound = length(extremeRts);
[xx, exponGain] = expGain(xx(:).');
if(nFound)
   raw(1:nFound) = extremeRts/exponGain;
end

% Arrays inside() and outside() contain the magnitudes of the most recently
% deflated roots inside and outside the unit circle. When I begin a search for
% a new root outside the unit circle, I use the values in outside() to determine
% the starting radius. Similarly for inside. In other words, I start at a radius
% where roots were recently found.
minOutside = 1 + 3/degree; % Start at least this far outside the unit circle.
saveLength = 5; % Number of previous radii to save.
inside  = 1 + zeros(1, saveLength); % To save recent root magnitudes inside the circle.
outside = minOutside + zeros(1, saveLength);
insideK  = 1; % Index to access inside()
outsideK = 1; % Index to access outside()
lastRadius = 0; % 1 means last root radius was > 1.   0 otherwise.

% Create the array of "random" arguments.

% I do not understand it, but it is MUCH better to have 2*pi and not pi
% in the definition of theta below. Using 2*pi seems wrong. It results in
% theta = 0, pi, pi/2, 3pi/2, pi/4, 5pi/4, 3pi/4, 7pi/4, ...
% The reason this seems wrong is that at unit magnitude, the third and fourth
% are complex conjugates of each other. So are the fifth and eighth and the
% sixth and seventh. When you do the third and fourth angle, you will deflate
% 2 roots (and their conjugates), that are extremely close to one another. This
% is highly non-random. The spectrum will be very large there. However, when 
% I use pi, on 10 mean-zero 10,000 degree polynomials, it failed every time. 
% When I used 2*pi on the same polynomials, the initial deflation found every 
% root without resorting to a final massive deflation. 

p = nextpow2(length(xx));
theta = 2*pi/(2^p)*bitReverse(p); % Quasi-random values between 0 and 2*pi

% theta(1)==0. If mean(xx) is significantly different from 0 then the spectrum is
% already large at 0 Hz. If you start a root search at theta(1) you will find a root
% near 1.0. When you deflate it, you boost the spectrum at 0 Hz which is already
% large. The resulting quotient can be enormous at 0 Hz. It will probably have
% a huge dynamic range and be ill-conditioned. Therefore, use theta(1) only if
% the mean is small.

% Consider two cases for x: 
%    x is random coefficient and non-negative or
%    x is random coefficient and mean 0.
% In the first case it is vitally important that you not deflate a root with 
% argument approximately 0. In the second case it is vitally important they you do
% deflate a root with argument approximately 0. Note theta(1:2)=[0,pi]. In the first
% case we want to skip theta(1) and start at theta(2).


% if(abs(mean(xx)) < .1) % This failed 5 times out of 5 tries at degree 100,000
if(abs(mean(normalizeJWF(xx))) < .1) % The previous normalizePow2 may have made mean(xx) large.
   thetaK = 1;                    % That is why I need normalizeJWF(xx) here.
else
   thetaK = 2;
end

% stopAt = 100;

while (nFound < degree)
%   if(abs(nFound-stopAt)<2)
%      pause;
%   end
%   if(rem(nFound,100) < 2)
%      disp(nFound)
%      plot(xx);
%      pause
%      plotdb(xx);
%      pause
%   end
   % Determine the radius to search for the next root.
   if(lastRadius == 1) % If last root was > 1, try to make next < 1
      r = median(inside);
   else
      r = max(minOutside, max(outside));
   end
   
   th = theta(thetaK); % Determine the angle to search for the next root.
   thetaK = thetaK + 1;
   
  % Don't use newton3 or laguerre3. Why? If the initial position is outside the unit
  % circle but not close to a root, they can converge to something far away.
  % Don't use newton2. I did once and got caught in the Zeno Effect.
   
   thisZ = laguerre2(xx, r*complex(cos(th),sin(th))); % Much better than Newton.
   absZ = abs(thisZ);
   if(absZ == Inf | absZ == NaN)
   %   fprintf('Failed to find a root at nFound=%g.\n',nFound);
   else
      lastRadius = (absZ > 1); % Was this root outside the circle? Yes or no.
      xx = deflate(xx, thisZ); % Deflate the root (and its conjugate if it's complex).
      
      if(absZ < 1) % Save this radius into inside() or outside() as appropriate.
         inside(insideK) = absZ;
         insideK = max(1, rem(insideK + 1, length(outside)+1)); % Increment insideK.
      else
         outside(outsideK) = absZ;
         outsideK = max(1,rem(outsideK + 1, length(outside)+1)); % Increment outsideK.
      end
          
      nFound = nFound + 1;
      raw(nFound) = thisZ;
      if(~isreal(thisZ)) % If it was complex, get its conjugate too.
         nFound = nFound + 1;
         raw(nFound) = conj(thisZ);
      end
   end
end

clear theta inside outside;

if(exponGain ~= 1)
   raw = exponGain*raw;
end

if(nargout > 1)
  % Use laguerre3 because of the stalking derivative root problem.
  [polished, err] = laguerre3(x, uhp(raw));
  [polished, err] = nonPolish(polished, err); % This also ensures they are uhp.
  [polished, err] = uniq(polished, err);
  fprintf('iterativeDeflation needs to find %g missing roots.\n', degree - hpdegree(polished));
  [polished,err,rc,outparms] = iterativeDeflation(x,polished,err,ID_NUM); % Deflate to find missing stragglers.
  [polished, err] = hp2fp(polished, err);
  if(rc < 0)
     fprintf('Unable to obtain a full set of accurate polished roots.\n');
  else
     fprintf('Return code=%g	Error=%g.\n',rc,outparms(2));
  end
end

if(NZeros) % If there initially were trailing zeros.
  raw = [raw; zeros(NZeros,1)];
  if(nargout >= 2)
     polished = [polished; zeros(NZeros,1)];
     if(nargout >= 3)
        err = [err; zeros(NZeros,1)];
     end
  end
end

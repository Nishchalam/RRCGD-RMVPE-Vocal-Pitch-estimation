% function [z, zErr, x, rc, parms] = extremeRoots(x, ID_NUM);
%
% PURPOSE
%   Find and deflate roots with extremely small or large magnitude. After doing that,
%   if there are only a few roots left, try to find all of them.
% 
% REQUIRRED INPUT ARGUMENTS
%   x - Any polynomial.
%   ID_NUM - id number for use in the history file
%
% REQUIRED OUTPUT
%   z - The extreme magnitude roots in the upper half-plane.
%
% OPTIONAL OUTPUT
%   zErr - Error estimates for z.
%   x - The polynomial after deflating these roots.
%   rc - See iterativeDeflation for the meaning.
%   parms - See iterativeDeflation for the meaning.
%   If there are too many roots to even try to find all of them, then
%   rc = parms(1) = parms(2) = -Inf.
%
% USES
%   bigRoots, flip, hp2fp, hpdegree, iterativeDeflation, lejaJWF, multroot, normalizePow2, 
%   polishMultRoot, polycmp

% AUTHOR
%   James W. Fox 2/24/04

function [z, zErr, x, rc, parms] = extremeRoots(x, ID_NUM);

if(nargin < 2)
   ID_NUM = 0;
end

sX = size(x);
if((sX(1)~=1 & sX(2)~=1) | ~isreal(x))
   error('$1 must be a 1 dimensional real array.');
end

xsave = x;
[zBig, zBigErr, x, rcBig] = bigRoots(x,1000);
x = normalizePow2(x);

% Here is, I think, an excellent idea for the future. zBig may have errors and the deflation
% process can introduce errors. Instead of flipping x, which has been deflated, write a new routine,
% smallRoots() and use the undeflated polynomial. Let z be the fourth roots of unity and let
% x=poly([1e-20*z 1e-20*z  1e10 1e10 1e10]). lroots(x) accurately finds the 3 large roots. However the 8
% small roots are found as 1.087770609288738e+026 with multiplicity 8. Not only is the multiplicity wrong,
% it is suppposed to be around 1e-20 and not 1e26! roots(x) is pretty good but multroot is a disaster.
% ==> Check is this is from  one of my modifications to multroot. If you use the undeflated poly,
% 1./multroot(flip(x(4:12))) is roughly ok but it does not find any multiple roots. However, if you use
% polishMultipleRoots and polish each with multiplicity 2, then you get outstanding results. 
% I.e., you need a program that can tell you the multiplicity of a root independent of multroot.
[zSmall, zSmallErr, x, rcSmall] = bigRoots(flip(x),1000);

zSmall = conj(1./zSmall); % They were upper half plane roots of the flipped polynomial.
z = [zSmall; zBig]; 
zErr = [zSmallErr; zBigErr];
x = normalizePow2(x);
rc = rcBig + rcSmall;
x = flip(x); % x was flippped by the previous line.
k=find(imag(z) >= 0);
if(~isempty(k))
   z = z(k);
   zErr = zErr(k);
end

parms = zeros(1,2);
parms(1) = length(hp2fp(z)); % May be increased below.
if(length(x) > 35) % Multroot would take too long.
   rc = -Inf;
   parms(1) = -Inf;
   parms(2) = -Inf;
else
   if(length(x) == 1) % Already found all the roots.
      parms(2) = polycmp(xsave, poly(lejaJWF(hp2fp(z))));
      return;
   end
   
   [z1d, z2d, parms1] = multroot(x,0,1e-9);
   [z1d, z1dErr] = polishMultRoot(x, z2d, parms1(1));
   [z1d, z1dErr] = uhp(z1d, z1dErr);
   [rts, rtsErr, rc1, parms1] = iterativeDeflation(x, z1d, z1dErr, ID_NUM);
   if(rc1 >= 0)
      x = 1;
      rc = rc + rc1;
      parms(1) = parms(1) + length(rts);
      z = [z; rts];
      zErr = [zErr; rtsErr];
      y = poly(lejaJWF(hp2fp(z)));
      if(length(xsave) == length(y))
         parms(2) = polycmp(xsave, y);
      end
   else
      rc = rc1;
   end
end

%multroot and pejroot can have problems. For example
%EDU» x=[1  1e-156  1e-136  0  -1e-80]; % Approximately [1 0 0 0 1e-80]
%EDU» multroot(x,0) % No, it is not a root of multiplicity 4
%ans =
%   -2.500000000000000e-157
%   -2.500000000000000e-157
%   -2.500000000000000e-157
%   -2.500000000000000e-157
% Using the flipped polynomial gives the answer we want.
% EDU» 1./multroot(flip(e),0)
%ans =
%                          0 +1.000000000000000e-020i
%     9.999999999999998e-021                         
%                          0 -1.000000000000000e-020i
%    -9.999999999999993e-021           
%
% Sometimes using pejroot, as in the section below, can improve accuracy. For example,
% x=poly(lejaJWF([1 1 1e11 1e11 1e13 1e13 1e15 1e15]));
% sort(lroots(x))
% ans =
%    9.999999999999998e-001
%    9.999999999999998e-001
%    1.000000000022189e+011
%    1.000000000022189e+011
%    9.999999999976357e+012
%    9.999999999976357e+012
%    1.000000000000450e+015
%    1.000000000000450e+015
% If you apply the section below you get almost perfect accuracy. 

% Also multroot(x) does an outstanding job. Perhaps you should use nroots and if the ratio
% of the largest to the smallest root is < 2e-16, then use multroot on the whole polynomial.

% However I have seen examples where the results are already quite accurate and pejroot
% significantly degrades the accuracy. One such polynmial is
% x = poly(lejaJWF([1e-100 1e-60 1e-30 1 1 1e20 1e20 1e40 1e40]));

% pejroot has problems. Sometimes the following section helps. Sometimes it hurts.
% Here is a big problem, either pejroot or laguerre3 can take a real root and give it an imaginary part.
% Then uhp() eliminates it. This is why many factorizations lost roots.

%if(rc >= 0 & hpdegree(z)==(length(xsave)-1)) % We found all the roots.
%   z = hp2fp(z);
%   z2d = z1d_z2d(z);
%   w = pejroot(xsave, z2d(:,1).', z2d(:,2)'); % Polish the roots.
%   w(:,1) = laguerre3(xsave, w(:,1)); % pejroot can turn 1.0, which is correct, to 1.0+1.18e-66
%   z = uhp(z2d_z1d(w));
%end
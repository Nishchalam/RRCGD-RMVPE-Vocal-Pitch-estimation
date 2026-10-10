% function [z, zErr, x, rc] = bigRoots(x,r);
%
% PURPOSE
%   Find and optionally remove roots with big magnitude.
% 
% REQUIRED INPUT ARGUMENTS
%   x - Any polynomial. x(1) must not be 0.
% 
% OPTIONAL INPUT ARGUMENTS
%   r - Default=100. If, nroots(x,r)==nroots(x,2*r) then return the
%       roots with magnitude > r. If not, let r=2*r and repeat the
%       previous sentence. It is suggested that r > 10, preferrably bigger.
%
% REQUIRED OUTPUT
%   z1d - The polished upper half-plane big magnitude roots. Must be a colum vector.
%
% OPTIONAL OUTPUT
%   z1dErr - Error estimates for z1d.
%   x - The polynomial after deflating z1d and its complex conjugates.
%   rc - The number of deflations
%
% USES:
%   flip, hp2fp, laguerre3, lejaJWF, nonPolish, normalizePow2, nroots, rootsUR

% AUTHOR
%   James W. Fox

% CHANGES
%   2/24/2004 Rewritten to use multroot and to use nroots more intensely.
%   2/24/2006 Ensure that we have exact complex conjugates.
%   3/3/2006 If the multroot error is small, skip using polishMultipleRoot

function [z, zErr, x, rc] = bigRoots(x,r);

xsave = x; % Polish against the original poly instead of x which is being deflated.

sX = size(x);
if(~x(1))
   error('$1(1) must not be 0.');
elseif((sX(1)~=1 & sX(2)~=1) | ~isreal(x))
   error('$1 must be a 1 dimensional real array.');
elseif(nargin==2)
   if(length(r)~=1 | ~isreal(r) | r <= 1)
      error('$2 must be a real scalar that is > 1.')
   end
end

degree = length(x) - 1;
nr = zeros(1,30);
for k=1:30
   % If I say nroots(x,1e10^k) then some user who wants it to fail will create
   % a polynomial with multiple roots at say 1e30. nroots can give the wrong answer
   
   % if you evaluate at exactly where there are multiple roots.
   nr(k) = nroots(x, mean(rand(1,10))*1e10^k); 
   if(nr(k) == degree)
      break;
   end
end
nr = nr(1:max(find(nr)));
d = diff(nr);
d = d(find(d));
z = [];
zErr = [];

if(length(nr) == 1)
   z1d = [];
   z1dErr = [];
   rc = 0;   
else
   rc = length(d);
   for k=rc:-1:1
      [z1d, z2d, parms] = multroot(x(1:(d(k)+1)),0,1e-9);
      
      if(max(z2d(:,2))==1 | parms(2)>1e-14) %3/3/2006 If there are no multiple roots or parms is large, use polishMultipleRoot
        % [z1d, z1dErr] = polishMultRoot(x, z2d, parms(1));
        % [z1d, z1dErr] = nonPolish(z1d, z1dErr, 0);
      
        % z=1e20 +1e20i;
        % x=poly([1e10 1e10 1e10 z conj(z) z conj(z) z conj(z)]);
        % format long e. It appears that z1d is perfect complex conjugates. However, the
        % imaginary parts differ in the last few bits. Use uhp and hp2fp to ensure they are perfect conjugates.
        z2d = uhp(z2d);
        z1d = z2d(:,1);
        z1dErr = z1d;
        sizeZ2d = size(z2d);
        for m = 1:sizeZ2d(1);
           % You would think it would be best to polish against the original poly instead of x which is being deflated.
           % However, amzaingly, tests show the opposite.
           % ==> Consider using pejroot instead of polishMultipleRoots <===
           [z1d(m), z1dErr(m)] = polishMultipleRoot(xsave,z2d(m,1), z2d(m,2)); % Often it is better to use x insted of xsave. Bizarre.
        end
      
        zWithErr = hp2fp_2d(z1d, z1dErr, z2d(:,2));
        z1d = zWithErr(:,1);
        z = [z; z1d];
        zErr = [zErr; zWithErr(:,2)];
      else
        z = [z; z1d];
        zErr = [zErr; parms(2)+zeros(length(z1d),1)];
      end
      
      % We are starting with the largest roots. Use backward deflation.
      numer = flip(normalizePow2(x));
      denom = flip(normalizePow2(real(poly(lejaJWF(z1d)))));
      if(~denom(1))
         error('bigRoots underflowed when unfactoring a set of roots.');
      end
%      x = (normalizePow2(bdeconv(numer,denom))); % This uses filter but is horribly slow.
      x = flip(normalizePow2(filter(1,denom,numer)));
      x = normalizePow2(x((end - length(numer) + length(denom)):end)); 
   end
end

% function z2d = hp2fp_2d(z, err, mul);
%
% PURPOSE:
%   Given a set of roots of a real polynomial (which may be only the upper half-plane) and root errors
%   and root multiplicities, create a nx2 array (roots,errors) with all roots in both half planes.
%
% REQUIRED INPUT
%  z - A set of roots which may only be upper half plane.
%  err - Error estimates for the roots
%  mul - Multiplicities of the roots
%
% REQUIRED OUTPUT
%  z2d - A 2d array of roots in both half planes and errors in column 2.

% AUTHOR
%  James W. Fox 2/28/2006

function z2d = hp2fp_2d(z, err, mul);

if(length(z) ~= length(err) | length(z) ~= length(mul))
   error('The 3 input args must have the same length.')
elseif(~isreal(err) | min(err) < 0)
   error('$2 must be a positive real array.')
elseif(~isreal(mul) | min(mul) <= 0 | ~min(mul == round(mul)))
   error('$3 must be a positive real integer array.')
end

% Restrict inputs to upper half plane.
k = find(imag(z) >= 0);
lenK = length(k);
if(lenK < length(z))
   z = z(k);
   err = err(k);
   mul = mul(k);
end

k1 = find(imag(z) == 0);
k2 = find(imag(z) > 0);
m = sum(mul(k1)) + 2*sum(mul(k2));
z2d = zeros(m, 2);
ix = 1;
for k = 1:length(z)
   for n = 1:mul(k)
      z2d(ix,1:2) = [z(k), err(k)];
      ix = ix + 1;
      if(imag(z(k)))
         z2d(ix,1:2) = [conj(z(k)), err(k)];
         ix = ix + 1;
      end
   end
end
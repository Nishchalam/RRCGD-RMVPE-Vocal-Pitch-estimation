% function z = clipSmallParts(z);
%
% PURPOSE
%  If the real or imaginary part of a root is very small relative to |z|, clip it to 0.
%
% INPUT ARGUMENTS
%  z - a two dimensional array with roots in column 1 and multiplicities in column 2.
%
% OUTPUT ARGUMENTS
%  z - the two dimensional array with clipped small parts.

function z = clipSmallParts(z);

s = size(z);
if(s(2) ~= 2)
   error('$1 must be a 2 dimensional array with roots in column 1 and multiplicities in column 2.');
end

for k = 1:s(1)
   zz = z(k,1);
   absZz = abs(zz);
   if(abs(imag(zz)) <= 1e-11*absZz)
      zz = real(zz);
   end
   if(abs(real(zz)) <= 1e-11*absZz)
      zz = imag(zz)*i;
   end
   z(k,1) = zz;
end
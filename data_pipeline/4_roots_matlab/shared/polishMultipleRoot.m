% function [w, werr] = polishMultipleRoot(x, z, n)
% 
% PURPOSE
%   Polish one approximate root with a known multiplicity.
%
% INPUT ARGUMENTS
%    x - A real coefficient polynomial.
%    z - A multiple root.
%    n - Default=1. The multiplicity of z.
%
% OUTPUT ARGUMENTS
%    w - The polished position.
%    werr - The maximum error in the root location.
%
% AUTHOR
%    James W. Fox 7/12/04

function [w, werr] = polishMultipleRoot(x, z, n)

if(nargin ~= 2 & nargin ~= 3)
   error('Supply 2 or 3 input arguments.');
elseif(~isreal(x))
   error('Currently the polynomial must have real coefficients.');
elseif(length(z) ~= 1 | ischar(z))
   error('$2 must be a real or complex constant, not an array or string.');
elseif(nargin == 2)
   n = 1;
elseif(length(n) ~= 1 | ~isreal(n) | n < 1)
   error('$3 must be a positive integer.')
end

n = round(n);
w = zeros(1,n);
NFound = 0;
for k = 1:n
   if(isreal(z))
      thisW = laguerre3(x,z);
      if(~isreal(thisW) & NFound == n-1)
         error('Found more than the specified number of roots.');
      end
      if(isreal(thisW))
         NFound = NFound + 1;
         w(NFound) = thisW;
      else
         NFound = NFound + 1;
         w(NFound) = thisW;
         NFound = NFound + 1;
         w(NFound) = conj(thisW);
      end
   
      if(NFound == n)
         break;
      end
   else
      thisW = laguerre3(x,z);
      NFound = NFound + 1;
      w(NFound) = thisW;
   end
   x = deflate(x,thisW);
end

meanW = mean(w);
werr = max(abs((w - meanW))); % This is too big, but I do not know what is a good value.
w = meanW;
% If laguerre returns a complex value and if you use deflate, you get rid of the root
% and its conjugate. It may work better to get rid of one and only one root at a time.
% However, if you deflate a single complex root, then x will have complex coefficients and
% you can not use laguerre3 because innerpart will complain.
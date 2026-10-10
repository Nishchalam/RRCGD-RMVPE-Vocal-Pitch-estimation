% function [z,j] = csort(z, method)
%
% PURPOSE:
%   Sort a complex array. Default = sort by abs(z).
%   Matlab's sort(z) is the same as csort(z) or csort(z,1).
%   If the roots are sorted by angle, then negative/positive purely real and
%   purely imaginary values are sorted by increasing magnitude.
%
% REQUIRED INPUT
%   z - any complex array
%
% OPTIONAL INPUT:
%   method - sort order. Default = 1
%          = 1 - sort by abs(z)
%          = 2 - sort by angle(z)
%          = 3 - sort by real(z)
%          = 4 - sort by imag(z)
%
% REQUIRED INPUT
%   z - the sorted values
%
% OPTIONAL INPUT:
%   j - indicies required to sort the original array.
%
% USES
%   Nothing.
%
% AUTHOR:
%   James W. Fox 11/18/2003. 

% CHANGES 
%   5/13/2004 add second output argument and sort purely imaginary values as well as purely real values.

function [z,j] = csort(z, method)

if(nargin == 1)
  method = 1;
else
  if(method~=1 & method~=2 & method~=3 & method~=4)
    error('$2 must equal 1, 2, 3, or 4.');
  end
end

if(method == 1)
  [junk,j] = sort(abs(z));
elseif(method == 2)
   ang = angle(z);
   ang(find(ang == -pi)) = pi; % Unfortunately, ang(k) can be -pi.
  [junk,j] = sort(ang);
elseif(method == 3)
  [junk,j] = sort(real(z));
elseif(method == 4)
  [junk,j] = sort(imag(z));
end

z = z(j);

if(method == 2) % Sort negative/positive reals and imaginaries.
   k = find(~imag(z)); % Indicies of purely real values
   if(~isempty(k))
      kk = find(z(k) < 0); % k(kk) are indicies of negative reals
      if(~isempty(kk))
         [junk,jj] = sort(z(k(kk)));
         z(k(kk)) = z(k(kk(jj)));
         j(k(kk)) = j(k(kk(jj)));
      end
   
      kk = find(z(k) >= 0); % k(kk) are indicies of positive reals
      if(~isempty(kk))
         [junk,jj] = sort(z(k(kk)));
         z(k(kk)) = z(k(kk(jj)));
         j(k(kk)) = j(k(kk(jj)));
      end
   end
   
   k = find(~real(z)); % Indicies of purely imaginary values
   if(~isempty(k))
      kk = find(imag(z(k)) < 0); % k(kk) are indicies of negative imaginaries.
      if(~isempty(kk))
         [junk,jj] = sort(z(k(kk)));
         z(k(kk)) = z(k(kk(jj)));
         j(k(kk)) = j(k(kk(jj)));
      end
   
      kk = find(imag(z(k)) >= 0); % k(kk) are indicies of positive imaginaries.
      if(~isempty(kk))
         [junk,jj] = sort(z(k(kk)));
         z(k(kk)) = z(k(kk(jj)));
         j(k(kk)) = j(k(kk(jj)));
      end
   end
end
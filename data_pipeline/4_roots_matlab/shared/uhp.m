% function [z, zerr] = uhp(z, zerr);
%
% Take a set of complex numbers and select out the values in the upper half plane.
% If there are also error estimates, select the corresponding values
%
% REQUIRED INPUT ARGUEMTS
%   z - any array of complex numbers
%
% OPTIONAL INPUT ARGUEMTS
%   zerr - an array of error estimates
%
% REQUIRED OUTPUT
%   z - The upper half plane roots
%
% OPTIONAL OUTPUT
%   zerr - The error estimates
%
% AUTHOR
%  James W. Fox

% CHANGES
%   2/28/06 Previously it assuemd that z was 1d. However, in bigroots I want to apply this to z2d which is 2d.
%           So change imag(z) to imag(z(:,1) and change z(k) to z(k,:).

function [z, zerr] = uhp(z, zerr);
k = find(imag(z(:,1))>=0);
if(~isempty(k))
   z = z(k,:);
   if(nargin == 2)
      zerr = zerr(k);
   end
end
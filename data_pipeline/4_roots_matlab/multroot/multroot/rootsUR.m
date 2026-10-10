% function z = rootsUR(x);
%
% PURPOSE
%   Matlab's roots() with underflow repair.
%
%   Use the eigenvalue method to find roots of x. 
%   Take corrective action if some roots underflow to 0.
%
% REQUIRED INPUT ARGUMENTS
%   x - Any polynomial.
%
% REQUIRED OUTPUT
%   z - The roots.
%
% USES:
%   nothing

% AUTHOR
%   James W. Fox 10/11/03

% This could be improved. Consider x=poly([1, 1, 1e32, 1e32]);
% Roots() underflows for the 2 small roots so rootsur factors x(3:5) and
% correctly finds that 1 is a double root. Now consider, x=poly([1, 1, 1e20, 1e20]);
% Unfortunately, roots() inaccurately finds 2.0 for one of the small roots and
% underflows for the other. Therefore, rootsur factors x(4:5) and gets 0.5. Not good.
% However if you factor x(3:5) then you correctly get 1 as a double root.
% Extremeroots() goes a long way to improving this.

function z = rootsUR(x);
   Lx = length(x);
   if(isempty(x))
      error('The empty polynomial does not have any roots.');
   end
   ixNZ = find(x ~= 0); % Indicies of non-zero values.
   if(isempty(ixNZ))
      error('Every complex number is a root of the zero polynomial.');
   end
   ixFNZ = ixNZ(1); % Index of the first non-zero.
   ixLNZ = ixNZ(end); % Index of the last non-zero.
   nZeros = Lx - ixLNZ; % Multiplicity of 0.0 as a root.
   if(nZeros)
      x = x(ixFNZ:ixLNZ); % Remove leading/trailing zeros so 0 is not a root.
   end

   z = roots(x);
   ixU = find(z == 0); % Indicies of underflows.
   nU = length(ixU);   % The number of underflows.
   while(nU)
      z(ixU) = roots(x(end-nU:end)); % Factor the last nU+1 coefficients.
      ixU = find(z == 0);
      nU = length(ixU);
   end
   
   if(nZeros)
      z = [zeros(nZeros,1); z]; 
   end

% function ixRelMin = ix2dmin(degree, x, y, z);
%
% PURPOSE:
%   Find indices of relative minima of the absolute values of a polynomial on 
%   a 2d region made up of 3 parallel slices. They may be 3 parallel circles. 
%   and if they are, they must not be complete circles.
%
% REQUIRED INPUT:
%   degree - the polynomial degree. Use to set the max length of the output.
%   x - The absolute polynomial values on one of the outer slices.
%   y - The absolute polynomial values on the inner slice.
%   z - The absolute polynomial values on the other outer slice.
%   You must have length(x) == length(y) == length(z) >= 3.
%
% OUTPUT:
%   ixRelMin - The indices where y is less than the surrounding 8 values.

% USES:
%   Nothing.
%
% AUTHOR:
%   James W. Fox 4/22/2001
%
%   Modified to run in constant memory by Ryan King, 12/9/2005
%   Modified to only save indices of relative minima instead of a much longer
%      array of 0's & 1's indicating the minima. Jim Fox, 12/10/2005
%   Modified to have 2 versions, one for Matlab 5 and one for Matlab 7. Jim Fox 12/10/05

function ixRelMin = ix2dmin(degree, x, y, z);

Lx = length(x);
if(length(y) ~= Lx | length(z) ~= Lx)
  error('All 3 arguments must have the same length.')
end

% Ensure that everything is a column.
x = x(:); % The output of is2dmin() will be used to construct
y = y(:); % an array of roots. The Matlab convention is that
z = z(:); % arrays of roots are column vectors.

ver = str2double(regexp(version,'^\d+','match','once')); % The major version of Matlab. Must be >= 5.

if(ver < 7)
   j = (2 : (Lx - 1)); 
   ixShiftLeft  = j-1; % Two long arrays which are used 3 times - cache them.
   ixShiftRight = j+1;

   thisy = y(j);
   temp =        (thisy < y(ixShiftLeft));
   temp = temp & (thisy < y(ixShiftRight));
   temp = temp & (thisy < x(ixShiftLeft));
   temp = temp & (thisy < x(j));
   temp = temp & (thisy < x(ixShiftRight));
   temp = temp & (thisy < z(ixShiftLeft));
   temp = temp & (thisy < z(j));
   temp = temp & (thisy < z(ixShiftRight));
   
   thisy = y(1);
   temp1 = thisy<x(1) & thisy<x(2) & thisy<y(2) & thisy<z(1) & thisy<z(2);
   
   thisy = y(end);
   temp2 = Lx*(thisy<x(end-1) & thisy<x(end) & thisy<y(end-1) & thisy<z(end-1) & thisy<z(end));
      
   temp = 1 + find(temp); % +1 because we started at y(2)
   Ltemp = length(temp);
   
   ixRelMin = zeros(Ltemp + length(find([temp1 temp2])),1);
   if(temp1)
      ixRelMin(1) = 1;
      ixRelMin(2:(Ltemp+1)) = temp;
   else
      ixRelMin(1:Ltemp) = temp;
   end
   if(temp2)
      ixRelMin(end) = temp2;
   end
else
   % Preallocate to avoid memory fragmentation.
   % There should only be degree output values, and the only time when 
   % you get this many is when the polynomial is x^degree-1.
   % However, we could get false positives so allocate a little bit more.
   ixRelMin = zeros(round(1.1*degree),1); % round() because it complains about fractions
   ix = 1; % Index into ixRelMin

   thisy = y(1); % The first term of y is a special case; there is no y(j-1)
   if(thisy < y(2) & ...
      thisy < x(1) & thisy < x(2) & ...
      thisy < z(1) & thisy < z(2))
   
      ixRelMin(1) = 1;
      ix = 2;
   end

   for j = 2:(Lx-1);
      thisy=y(j);
      if(thisy < y(j-1) &                thisy < y(j+1) & ...
         thisy < x(j-1) & thisy < x(j) & thisy < x(j+1) & ...
         thisy < z(j-1) & thisy < z(j) & thisy < z(j+1))
         % && is better than & for tests on multiple scalars.
      
         ixRelMin(ix) = j;
         ix = ix + 1;
      end;
   end

   thisy = y(Lx); % The last term of y is a special case; there is no y(j+1)
   if(thisy < y(Lx-1) & ...
      thisy < x(Lx) & thisy < x(Lx-1) & ...
      thisy < z(Lx) & thisy < z(Lx-1))
   
      ixRelMin(ix) = Lx;
      ix = ix + 1;
   end
   ixRelMin = ixRelMin(1:ix-1); % Resize to the values actually found
end

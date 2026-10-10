% function  plotroots(z1, z2, z3)
%
% PURPOSE
%   Plot 1, 2, or 3 sets of complex numbers in polar coordinates.
%   If there is more than 1 array, they can have different lengths.
%
%   Only values in the upper half plane will be plotted because it
%   is assumed the numbers are roots of real polynomials and thus
%   the values in the lower half plane are symmetric.
%
%   The first set uses blue x's. The second uses red o's. 
%   The third uses green diamonds.
%
% AUTHOR
%  James W. Fox

function  plotroots(z1, z2, z3)

LOGRADIUS = 0; % Non-zero to use log(radius) instead of radius for y axis.
if(LOGRADIUS)
   UC = 0; % Ordinate of the x-axis which corresponds to the unit circle.
else
   UC = 1;
end

z1 = z1(find(imag(z1) >= 0));
if(LOGRADIUS)
   r1 = log(abs(z1));
else
   r1 = abs(z1);
end
a1 = angle(z1) * 180/pi; % Radians to degrees
rmin = min(r1);
rmax = max(r1);
if(nargin == 1)
   plot(a1,r1,'bx', [0,180],[UC,UC],'r');
else
	z2 = z2(find(imag(z2) >= 0));
   if(LOGRADIUS)
      r2 = log(abs(z2));
   else
      r2 = abs(z2);
   end
   a2 = angle(z2) * 180/pi;
   rmin = min(rmin, min(r2));
   rmax = max(rmax, max(r2));
   if(nargin == 2)
      plot(a1,r1,'bx', a2,r2,'ro', [0,180],[UC,UC],'r');
   else
		z3 = z3(find(imag(z3) >= 0));
      if(LOGRADIUS)
         r3 = log(abs(z3));
      else
         r3 = abs(z3);
      end
		a3 = angle(z3) * 180/pi;
	   rmin = min(rmin, min(r3));
   	rmax = max(rmax, max(r3));
	   plot(a1,r1,'bx', a2,r2,'ro', a3,r3,'gd', [0,180],[UC,UC],'r');
   end
end   
title('Plot of Roots in Polar Coordinates')
xlabel('Angle in degrees');
if(LOGRADIUS)
   ylabel('log(Radius)');
else
   ylabel('Radius');
end
   % If rmax-rmin=4E-16 you get an error message from plot() and axis().
% If I reset rmax and rmin, I still get an error message from axis()
% but the axis call works none the less. However, if I do nothing,
% I get an error messge but Matlab does something intelligent anyway.
axis([0, 180, rmin, rmax]); % Ensures we go from 0 to 180 degrees.

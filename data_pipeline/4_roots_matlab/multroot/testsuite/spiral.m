function [p,z] = spiral
%
%  test polynomial used by James W. Fox. Roots are on a spiral plus a double root.
%
   n=60;
   z=(1.00001.^(0:(n-1))) .* exp(i*(2*pi/n)*(0:n-1));
   p=poly(leja(z)); % A poly with roots on a spiral
   p=conv(p,[1 -1]); % Make 1 a double root.
   z = [z.',ones(n,1)];
   z(1,2) = 2; % The first root is 1 and it is a double root.
   z = displayRoots(z,2); % This also clips small parts and sorts

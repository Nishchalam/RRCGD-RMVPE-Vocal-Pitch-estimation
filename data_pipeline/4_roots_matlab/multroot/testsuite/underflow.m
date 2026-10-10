function [p,z] = underflow
%
%  test polynomial used by James W. Fox. A cubic polynomial,
%  p = 1.0000e+000 -1.0000e+032  2.0000e+032 -1.0000e+032.
%  Matlab and the original multroot underflow for the two smaller roots.
%  The newer multroot correctly finds that 1 is a double root.
   z = [1e32, 1, 1];
   p = poly(z);
   z = [z.',ones(3,1)];
   z = displayRoots(z,2); % This also clips small parts and sorts

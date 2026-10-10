% function tnroots(x, radii,z)
%
% PURPOSE
%  Test nroots
%
% AUTHOR
%  James W. Fox

function tnroots(x, radii,z)
% radii - an array of radii to check
% z - the roots of x. Used to get the correct answer
absZ =abs(z);
for k=1:length(radii)
   r = radii(k);
   correct = nels(ix(absZ<r));
   disp(sprintf('\n====================='));
   disp(sprintf('\nAt r=%g the correct answer is       %g',r,correct));
   disp(sprintf('Computer answer from winding number = %d',nroots(x,r))); 
   disp(' ');
   disp(sprintf('Computer answer from nrootsB =        %d',nrootsB(x,r))); 
   pause;
end
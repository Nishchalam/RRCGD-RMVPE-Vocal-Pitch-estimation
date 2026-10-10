% function rts = quadraticRoots(coeff)
%
% PURPOSE
%   Find the roots of a quadratic equation.

% USES
%   normalizePow2
%
% AUTHOR
%   James W. Fox 10/6/2003

function rts = quadraticRoots(coeff)
if(length(coeff) ~= 3)
   error('Must have length($1) == 4');
elseif(coeff(1) == 0)
   error('Must have $1(1) ~= 0.');
end

% I once had a poly with coefficients: 
% 7.466108948025740e-301 3.101300322905030e-291 7.466108948025739e-251
% The discriminant underflowed to 0. Normalizing cured this.
coeff = normalizePow2(coeff);

a = coeff(1);
b = coeff(2);
c = coeff(3);
discriminant = sqrt(b*b - 4*a*c);
q = -(b + sign(b)*(discriminant))/2;
% Note: if b==0 then q==0.
if(abs(q) > 1e-13) % Do not divide by zero.
   rts = [q/a; c/q]; % The roots of the quadratic computed stably.
else
   rts = [-b + discriminant; -b - discriminant]/(2*a);
end

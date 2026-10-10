% function a = angl(z)
%
% PURPOSE
%   Return complex arguments in the range -pi < angle <= pi.
%   Note, Matlab returns them in the range -pi <= angle <= pi.
%
% INPUT ARGUMENTS:
%   z - Any complex array.
%
% Matlab probably uses atan2() to compute angle() and "help atan2()" says 
% the range of the answer is -pi <= and <= ans.
%
% Note the following from a Matlab fragment:
%>> x=[-1, i];
%>> angle(x) % I agree with this answer.
%	 ans = 3.1416    1.5708
%>> angle(conj(x)) % I DISAGREE with this answer.
%	 ans = -3.1416   -1.5708
% conj(-1)==-1 and yet angle(x) returns pi for this value while
% angle(conj(x)) returns -pi. If you had said x=[-1] then
% angle(conj(x)) correctly returns pi. It is only if x has both
% real and complex values that this anomaly occurs.
%
% I suspect that what is going on is the conj() function looks to see
% if any values have non-zero imaginary parts. If so, it changes the
% sign bit for all imaginary parts. If there was initially some real
% value, like -1, it is implicitly -1 + 0*i. Changing the sign bit
% gives -1 -0*i. angle(-1+0*i) returns pi while angle(-1-0*i) returns -pi.
% Typing angle(-1-0*i) does NOT return -pi but I suspect that this is 
% in essence what is going on above.

% WRITTEN BY:
%   James W. Fox 11/10/2001

function a = angl(z)

a = angle(z);
k = find(a == -pi);
if(length(k))
   a(k) = pi;
end
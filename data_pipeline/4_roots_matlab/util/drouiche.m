% function b=drouiche(n)
%
% PURPOSE
%   Compute the bit reversal of 0:(2^n-1) by Drouiche's method.
%   Karim Drouiche, "A New Efficient Computational Algorithm for Bit-Reversal Mapping,"
%   IEEE Transactions on Signal Processing, V 49, N 1, Jan 2001.

% AUTHOR
%  James W. Fox
%
% COMMENT
%   There are two differences with his paper. First, he initialized the first 3 elements. It
%   is actually only necessary to ensure the first element is zero. Also, in his indexing
%   the first element has index 0 while Matlab requires that the first element have index 1.
%   Also, his C code has bugs.

function b=drouiche(n)
if(n > 25)
   error('$1 must be < 26'); % This is true on my computer. JWF.
end
b = zeros(1,2^n); % Allocate output array
up = 1;     % Will be increased by succesive powers of 2.
down = 2^n; % Will be decreased by succesive powers of 2.
for k=0:(n-1)
   lastup = up;
   up     = 2*up;   % Equivalently, up = bitshift(up,1)
   down   = down/2; % Equivalently, down = bitshift(down,-1)
   b(up:-1:lastup+1) = (b(lastup) + down) - b(1:lastup); 
end
% You can use bitxor instead of "-" in the last line, but it is MUCH slower.

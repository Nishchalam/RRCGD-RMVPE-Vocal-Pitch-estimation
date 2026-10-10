% function b = bitReverse(n);
%
% Compute the bit reversal of 0 : (-1 + 2^n)
% bitReverse(n)/2^n is the first 2^n terms of the van der Corput sequence.
%
% This is useful for Monte Carlo methods, numerical quadrature,
% and for computing the fast Fourier transform.
%
% There are two terms used with random numbers - pseudo-random and quasi-random.
%
% A sequence of pseudo-random numbers is computed by some formula and
% is thus not truly random. However, statistically speaking, it is
% indistinguishable from a true sequence of random numbers. It passes
% statistical tests that true random numbers would pass.
%
% Quasi-random numbers are numbers that are highly uniformly distributed.
% They are desirable in Monte Carlo methods where true random sampling
% is less important than nearly uniform sampling. The elements of 
% bitReverse(p), in order, fill out 0 : (-1+2^n) in a nearly uniformly
% distributed manner. I.e., if you stop after k elements, they are
% nearly uniformly distributed. If k is a power of 2, they are exactly
% uniformly distributed. It has been proved that using quasi-random numbers 
% produces less error in the numerical quadrature of Monte Carlo than 
% pseudo-random numbers.
%
% The van der Corput sequence is a quasi-random sequence between 0 and 1.
% There are other sequences which, on average, are even more uniformly distributed.
% However, as the code below shows, the van der Corput sequence is extraordinarily
% quick to generate. The only draw back of the code below is that you must 
% allocate an array of length 2^n. If you are going to do a Monte Carlo simulation
% with a billion terms, you probably do not want an array with a billion terms.
%
% REQUIRED INPUT ARGUMENTS
% ------------------------
% n - Any non-negative integer. 
%
% OUTPUT ARGUMENTS
% ----------------
% b - The bit reversal of [0 : (2^n - 1)].
%
% USES
% ----
% Nothing
%
% AUTHOR
% ------
% James W. Fox 3/31/2003.

% function drouiche(n) does the same thing, with the same median times but the average
% time is 7% faster than bitReverse(). I doubt if the same is true for C.

% This algorithm was discovered while analyzing Drouiche's paper on bit-reversal.
% It has the same complexity as Drouiche but it is easier to understand. Also, if it
% is converted to C, it may be faster if you use C's bitwise "or" instead of "+."

function b = bitReverse(n);
if(n > 25)
   error('$1 must be < 26'); % At least this is true on my computer. JWF.
end
b = zeros(1, 2^n); % Allocate the output array.
for k = 0:(n-1)
   b((2^k+1):2^(k+1)) = b(1:2^k) +  2^(n-k-1);
end

% Here are two previous versions. This is better for translation to C.

%n = 2^p; % n = bitshift(1,p);
%topDown  = n; % Successively takes on the values n/2, n/4, n/8, ...
%bottomUp = 1; % Successively takes on the values 1, 2, 4, 8, 16,...
%r = zeros(1, n); % Allocate space for the answer.

%for k = 1:p
                                %%    j = [1:bottomUp];
%    topDown = topDown/2;       %% topDown = bitshift(topDown,-1);
%    r(bottomUp+1:2*bottomUp) = topDown + r(1:bottomUp); %% r(bottomUp+j) = bitor(topDown,r(j));
%    bottomUp = 2*bottomUp;     %% bottomUp = bitshift(bottomUp,1);
%end

% The double commented out code, %%, also works but is 5 times slower. 
% However, note the use of bitshift() and bitor() which may be faster
% in C than addition and division.
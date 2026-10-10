% function [y,err] = innerpart(x, zin)
%
% PURPOSE
%   Compute the inner part of a polynomial evaluation at an array of locations.
%
%   If the mex file innerpart.c is compiled, the following will never be used
%   and this is good because it is very slow if you ask for error estimates.
%
% INPUT ARGUMENTS
%   x - the polynomial
%   zin - an array of points to evaluate
%
% OUTPUT ARGUMENTS
%   y - The inner part of the polynomial evaluations at zin.
%   err - Error estimates for each evaluation.

% Suppose |z|>1. Then mathematically, innerpart(x,z)==innerpart(flip(x),1/z).
% That will be true for the code below. It is not true for innerpart.dll.
% The reason is that when I compute 1/z in C, I get a different answer than
% Matlab. It is common to differ by 3 in the 16th digit.

% Do not use this. Instead compile the mex file innerpart.c. Err is set to 1e-7
% below. That is OK for random coefficient, but not for a general polynomial.
% The mex file correctly computes error estimates.

% AUTHOR:
%   James W. Fox 4/22/2002

function [y,err] = innerpart(x, zin)

y   = zeros(size(zin));
err = y;
if(max(abs(zin)) > 1)
   flippedX = flipud(x(:));
end

absZ = abs(zin);

k = find(absZ <= 1);
[y(k),err(k)] = hornerWithError(x,zin(k));

k = find(absZ > 1);
if(~isempty(k))
   [y(k),err(k)] = hornerWithError(flippedX,1./zin(k));
end

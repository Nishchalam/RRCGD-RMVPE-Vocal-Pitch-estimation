% function s = lejaJWF(z)
%
% PURPOSE:
%   Leja sort an array of complex numbers.
%
% REQUIRED INPUT:
%   z = Any 1D array of complex numbers
%
% OUTPUT:
%   s = The sorted values.
%
% DEFINITION:
%   The first value in a Leja sorted array is the complex value with the
%   largest absolute value. The second value is the value furthest from the
%   first value. Suppose n values have been sorted. Take each value that
%   has yet to be sorted and compute the product of the distances to the
%   sorted values. The value with the largest product becomes the next
%   sorted value. In other words, at each step you find the value that,
%   collectively, is furthest from all the previously sorted values.
%
% ALGORITHM:
%   Suppose there are n already sorted values. For each unsorted value
%   we must compute the product of the distances to the sorted values - a
%   product with n terms. However, on the previous iteration we computed
%   the product of the distances to the first n-1 sorted values. There
%   is no reason to repeat this product of n-1 values. Just save it and
%   reuse it. If we save and reuse products, we can compute the new 
%   product with a single multiplication instead of n multiplications.
%
% COMMENT 1: 
%   Leja ordering a set of polynomial roots before applying poly can
%   dramatically improve the results. For example, let 
%   z=exp(2*pi*i*(1:100)/100);
%   Then max(abs(poly(z))) is 8.9E7! However, poly(lejaJWF(z)) has 1 and
%   -1+1e-14i for the first and last value and the maximum absolute value
%   of the other values is 6e-15 (ideally they would be 0).
%
% COMMENT 2:
%   There is no such thing as THE Leja order. There can be more than one
%   correct answer. Leja talked about maximizing the products of distances
%   but he seems to not have considered what to do in the case of a tie. 
%
%   Example 1: lejaJWF(1:5) can be either [5 1 3 2 4] or [5 1 3 4 2].
%   Note, the product of the distances from 2 to [5 1 3] is 3*1*1 while the
%   product of the distances from 4 to [5 1 3] is 1*3*1. 4 and 2 tie.
%
%   Example 2: Suppose the value with largest absolute value is complex 
%   and suppose its complex conjugate is also in z. Then either one can be
%   the first value. They tie.
%
%   Example 3: Suppose the largest value is real and suppose the value 
%   furthest from it is complex. Suppose its conjugate is also in z. Then 
%   either one of these 2 values can be the second value. They tie.
%
%   Example 4: I have seen several instances where there were 2 different
%   equally valid Leja sorts. However, in every case, when values differed,
%   one value was the complex conjugate of the other.
%
% COMMENT 3: 
%    I am not sure if anyone has addressed the problem of how to sort 
%    multiple roots. For example, suppose you want to sort [1,1,2,2,3,3].
%    The first 3 values must be [3,1,2]. However this leaves unsorted 
%    [1,2,3] and in every case the product of the distances to the sorted
%    values is 0. We can not pick the maximum product, as in the definition
%    of Leja sorts, because all products are the same - they are 0.
%
%    This function adopts the following strategy. It first does a Leja sort
%    of all the unique values. Then, if there are multiple roots, it calls
%    itself recursively to sort the multiple roots. If n is the maximal
%    multiplicity, it will go to a recursion depth of n-1. For example,
%    this function will sort [1,1,1,2,2,2,3,3,3] as [3,1,2,3,1,2,3,1,2]
%    because the Leja sort of [1,2,3] is [3,1,2].
%
%    There is another Matlab code on the web to do Leja sorting. It has an
%    easily correctible bug, but it does not do anything to accomodate
%    multiple roots and this can be disasterous. For example let z be the
%    100'th roots of unity as in Comment 1. If you use this Leja to sort
%    [z,z] and then apply poly, the maximum absolute value is 8.5E7. Why?
%    The first 100 values in the sorted array are the 100'th roots of
%    unity in Leja order. However, the last 100 values are z exactly. By
%    contrast, for the values that should be 0 in poly(lejaJWF([z,z])) the
%    maximum absolute value is 1e-14.

% USES:
%   Nothing
%
% AUTHOR:
%   James W. Fox 8/19/2005

function s = lejaJWF(z)

sz = size(z);
if(prod(sz) <= 1) % Emptyset or singleton.
   s = z;
   return;
elseif(sz(1)~=1 & sz(2)~=1) % Mathworks might want to allow 2D arrays which
   error('$1 must be 1xn or nx1'); % are sorted by columns.
end

s = zeros(size(z)); % The output array will be the same shape as z.
indicies = s; % z(indicies) is in Leja order.
[junk,ixMax] = max(abs(z)); % If a tie, max returns the index of the first.
nextValue = z(ixMax); % The largest absolute value.
nFound = 1;
indicies(1) = ixMax;
s(1) = nextValue;  % The first sorted value.
products = abs(z - nextValue); % The product of distances to sorted values.
products(ixMax) = 0; % Make sure this index will never again hold the max.

for k=1:length(z)-1
   [maxV,ixMax] = max(products);
   if(~maxV) % This happens if only multiple roots are left to be sorted.
      ix = 1:length(indicies);
      ix(indicies(1:nFound)) = 0;  % Set indicies of sorted values to 0.
      ix = find(ix); % Indicies of unsorted values
      s(k+1:end) = lejaJWF(z(ix)); % Recursively sort the unsorted values.
      return
   end
   nextValue = z(ixMax);
   nFound = nFound + 1;
   indicies(nFound) = ixMax;
   s(nFound) = nextValue; % The next sorted value.

   % There should be no reason to do lejaJWF(rand(1,500)) but if you do,
   % there can be a problem. The distance between any two points is < 1. 
   % When you compute the products of the distances to already sorted 
   % values, this will be the product of SEVERAL values < 1. This may 
   % underflow. The following normalization will take care of this case. 
   % However, it will not work if you Leja sort 1e-5*rand(1,500). 
   % I.e., the numbers are much smaller to start with. You would need to do
   % much more frequent normalizations. You might therefore think it wise 
   % to normalize z at the beginning. That will indeed take care of 
   % 1e-5*rand(1,500) but then I have seen lejaJWF(rand(1,500)) overflow! 
   % The following will not take care of all underflow/overflow issues, but
   % it helps. It would be safer to do more frequent normalizations but 
   % that would slow the program down and in normal situations you will 
   % never need any normalizations.

   % Warning, do not move the next 3 lines after "products(ixMax)=0;"
   % Why? Suppose you have a polynomial of exactly degree 100. Then when
   % k=100 you will set the last value of products to zero. Thus
   % max(products)==0 and you will then divide by 0.
   if(~rem(k,100))
      products = products/max(products);
   end

   products = products .* abs(z - nextValue); % Update products.
   products(ixMax) = 0; % Make sure this index will never again be the max.
end

% function ans = permuteVCD(x);
%
% PURPOSE:
%   Permute the values of x by bit reversing its indicies.
%   "VDC" stands for Van der Corput who wrote about the
%   usefullness of bit reversal.
%
% INPUT:
%   x = Any 1D array.
%
% OUTPUT:
%   ans = The permuted array.

% USES:
%   bitReverse
%
% AUTHOR:
%   James W. Fox 4/26/2001

function ans = permuteVCD(x);

Lx = length(x);
if(Lx)
  j = bitReverse(ceil(log2(Lx)));
  j = 1 + j(find(j < Lx)); % "1+" since bitrevf() starts at 0.
  ans = x(j);
else
  ans = [];
end

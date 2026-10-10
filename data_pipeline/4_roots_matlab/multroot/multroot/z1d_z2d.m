% function [z2d] = z2d_1d(z1d)
%
% PURPOSE
%   Convert a 1d array of roots to a 2d array. The 2d array has unique roots in 
%   column 1 and their multiplicities in column 2.
%
% INPUT ARGUMENTS
%   z1d - a 1d array of roots
%
% REQUIRED OUTPUT
%   z2d - a 2d array of roots where column 2 is the multiplicity.

% AUTHOR
%   James W. Fox 06/22/04.

function [z2d] = z1d_z2d(z1d)

if(min(size(z1d)) ~= 1)
   error('%1 must be a 1 dimensional array.');
end

u = unique(z1d);
z2d = zeros(length(u),2);
z2d(:,1) = u;
for k = 1:length(u)
   z2d(k,2) = length(find(z1d == u(k)));
end

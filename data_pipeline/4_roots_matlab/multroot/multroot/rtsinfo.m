function [z,bkerr,cond] = rtsinfo(p)

% roots, backward error and condition number using "roots"

n = length(p)-1;

% JWF
% rootsUR can often find roots that incorrectly underflow to 0.           
% z = [roots(p),ones(n,1)];
z = [rootsUR(p),ones(n,1)]; % roots() with underflow repair.

p = p/p(1);
cond = spcond(z);

% JWF
% The next line is wrong because z is 2d, not 1d. If you use the old code this causes
% multroot([1, -5]) to compute horribly wrong error estimates. Why? The old code expects 
% that the sort will be applied to each column separately - the first column containing
% unique roots and the second column containing multiplicities. This is true provided there
% are at least 2 rows. However, if p=[1,-5] then z=[5,1] and Matlab does not sort this by
% columns, it sorts it as [1,5] because it sorts row vectors differently than it sorts
% matrices with more than 1 row.
% Also, multiplying together the roots of a polynomial can be much more accurate if you
% first leja sort them. Also, the new code is easier to understand - one line to unfactor
% the roots instead of a for loop.
%
%[yy,jj]=sort(abs(z)); y = z(jj);
% The correct line would be:
%[yy,jj] = sort(z(:,1)); y = z(jj,1);
%
%ff = [1]; w = ones(1,n+1);
%for k = 1:n
%    ff = conv(ff,poly(y(k)));
%    if abs(p(k+1)) > 1 
%        w(k+1) = 1/abs(p(k+1));
%    end
%end
ff = poly(lejaJWF(z(:,1))); % Unfactor all the roots.
w = ones(1,n+1); % Find the weights.
k = find(abs(p) > 1);
if(~isempty(k))
   w(k) = 1./abs(p(k));
end

bkerr = norm( (ff-p).*w, inf);

function cond = ...
                 spcond(z)
% SPCond calculates the struture preserving condition number
% of roots z(:,1) with multiplicity structure z(:,2)
% 
%   Input:  z --- mx2 matrix, z(:,1) is the root vector
%                             z(:,2) is the multiplicity structure
%
%  Output:  cond -- the structure preserving condition number
%
% Syntax:  >> cond = spcond(z)
%
   m = size(z,1);   % number of distinct roots
   n = sum(z(:,2)); % degree of the polynomial
   
% JWF
% Multiplying together the roots of a polynomial can be much more accurate if you
% first leja sort them. Also, the new code is easier to understand - 2 lines versus 8 lines.
% Also, it is possible for z(:,1) to contain 0.0 so sort(1./abs(z(:,1))); could 
% produce a divide by zero error. See an example of this below. The new 
% code will produce the same sort but without the possibility of division by zero.
% Here is an example of how to generate the error using the old version of multroot:
% EDU» multroot([1 -2 1 0])
% Warning: Divide by zero.
%> In c:\!mydocs\!matlab\multroot\multrootOrig\spcond.m at line 17
%  In c:\!mydocs\!matlab\multroot\multrootOrig\multroot.m at line 82
% The new code fixes this problem.
%
% sort
%  [yy,jj] = sort(1./abs(z(:,1)));   y = z(jj,1); l = z(jj,2)';
%   
%  ff = [1];  ll = l - ones(1,m);    % form the base polynomial
%  mm = max(ll); lll = ll;
%  for kk = 1:mm
%      id = find(lll>0);
%      ff = conv(ff,poly(y(id)));
%      lll = lll - 1;
%  end;
   [junk,jj] = sort(abs(z(:,1))); jj=flipud(jj); y=z(jj,1); l=z(jj,2);
   ff = poly(lejaJWF(z2d_z1d([y,(l-1)])));

   for j = 1:m               % remaining polynomial multiplication
       if m > 1
% JWF
% Multiplying together the roots of a polynomial can be much more accurate if you
% first leja sort them.
%          id = [1:(j-1),(j+1):m]; g = poly(y(id));
           id = [1:(j-1),(j+1):m]; g = poly(lejaJWF(y(id)));
       else
            g = [1];
       end;
       g = conv(-l(j)*g,ff);  J(:,j) = conj(g');
   end
   g = conv(-g/l(m), [1,-y(m)]);  % evalulate the coef. operator
   c = conj(g(2:n+1)');
  
   for j = 1:n                    % scale J
       if abs(c(j)) >= 1.0 
           s = abs(1.0/c(j));
           J(j,:) = J(j,:)*s;
       end
   end;
   
   [s,x] = zminsv(J,1.0e-16);       % finding smallest singular value
   
   cond = 1/s;
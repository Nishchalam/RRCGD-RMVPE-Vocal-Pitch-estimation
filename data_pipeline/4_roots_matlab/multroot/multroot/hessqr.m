
function [B,t] = hessqr(A)
%
%  assume A is a Hessenberg matrix
%  output  B -- upper triangular
%          c -- rotation used
%
   [m,n] = size(A);
   if m < n, return; end;
   
   B = A;
   
   for j = 1:n
       if j < m
           d = sqrt(B(j,j)^2+B(j+1,j)^2);
% JWF. Suppose x = poly(leja([1e-100 1e-60 1e-30 1 1 1e20 1e20 1e40 1e40]));
% Then the original multroot(x) complains that the following line is division by 0 but
% it still gave an acceptable answer. My revised version of multroot dies in pejroot. This
% is the only example I know where my revised version does worse than the original.
           c = B(j,j)/d; s = B(j+1,j)/d;
           T = [c,s;-s,c];
           B(j:j+1,j:n) = T*B(j:j+1,j:n);
           t(1,j) = c; t(2,j) = s;
       end;
   end;
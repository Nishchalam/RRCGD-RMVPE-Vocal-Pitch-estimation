% function [zout, zErr, nIter] = newton2(x, zin, maxNIter);
%
% This is the original polish.m. If |z|>1 it applies a correction 
% to a point outside the unit circle. It is essentially classical Newton.
% It allows the Zeno effect and stalking derivative roots.
%
% PURPOSE:
%   Newton's method optimized for polynomials. In particular, polynomial
%   evaluations never overflow, no matter how large the root or how high
%   the degree. If it does not converge in maxNIter steps, Inf is returned.
%
% REQUIRED INPUT:
%   x - polynomial coefficents. x(1) is the constant term.
%   zin - an array of starting points.
%
% OPTIONAL INPUT:
%   maxNIter - Default = 250. Maximum number of iterations.
%             If you do not mind returning Inf because the caller
%             will restart at a new position, then set maxNIter
%             small. If you have found all the roots and you wish
%             to polish their positions, set maxNIter large.
%             I have seen simple examples where 200 is too small.
%
%             maxNIter provides one kind of step size control which
%             will keep you from wandering too long in the fractal zone.
%             You may wish to implement additional step size control to
%             keep from correcting far outside the unit circle because
%             it will take a long time to come back. An intelligent 
%             default for maxNIter probably depends on the degree, 
%             but I do not know what it is.
%
% REQUIRED OUTPUT:
%   zout - A column vector of polished roots. Inf at locations where it
%          failed to polish in maxNIter steps.
%
% OPTIONAL OUTPUT:
%   zErr - An error estimate. Currently, this is the final polishing correction.
%   nIter - The number of iterations used for each root.
%
% ALGORITHM:
%   Let f be the polynomial and df be its derivative. Each iteration step does:
%
%   if(abs(z) <= 1)
%     z = z - polyval(f,z)/polyval(df,z)
%   else
%     z = z - z*polyval(flip(f),1/z)/polyval(flip(df),1/z)
%   end
%
% JUSTIFICATION:
%   Suppose f is an n degree polynomial and df is its derivative. 
%   The usual Newton iteration step is:
%	 	z = z - polyval(f,z)/polyval(df,z)
%   Note polyval(f,z) = z^n * polyval(flip(f), 1/z) and
%   polyval(df,z) = z^(n-1) * polyval(flip(df), 1/z)
%   Thus Newton can be rewritten as 
%	 	z = z - z*polyval(flip(f),1/z)/polyval(flip(df),1/z)
%   What is the point? If |z| >> 1, then polyval(f,z) may overflow 
%   while polyval(flip(f), 1/z) will not overflow. In fact it is 
%   bounded above by the sum of the absolute values of the coefficients.
%   Thus this routine will never return Inf because the polynomial
%   evaluation overflowed. (However, it can return Inf if it 
%   failed to converge in the allowed number of steps.)

% USES:
%   innerpart, normalizepow2, quadraticRoots
%
% WRITTEN BY:
%   James W. Fox, 4-03-2001

% Strange things can happen. Let w=poly(1:20) then polyval(flip(w),.05)==0.
% polyval(flip(w),.04999990026458411) == 8.970602038971265e-014
% polyval(flip(w),.04999990026458412) == 0
% polyval(flip(w),.04999990026458413) == 7.505107646466058e-014
% If you use horner() with output error estimate on the evaluation, in each case
% the error estimate is 1e-12 and the polynomial values are 8.97e-14, 0, and 7.5e-14 so 
% this criterion would have told you to stop at any of these. However, 
% 1/.04999990026458412==20.00003989424593 is not that great of a root estimate.

function [zout, zErr, nIter] = newton2(x, zin, maxNIter);

DIAGNOSTICS = 2; % 0=none, 1=minimal, 2=full diagnostic messages

if(nargin < 2)
   error('Supply at least 2 arguments.');
elseif(~length(x)) 
   error('Polynomial is empty.');
elseif(~any(x))
   error('Polynomial is entirely zero.');
end

if(x(1) == 0) % Remove leading zeros for speed and for linear polynomials.
   x = x(ixFirstNZ(x):end);
end

x = normalizepow2(x(:).'); % Matlab convention is that a polynomial is a row vector.
degree = length(x) - 1;
zout = zeros(size(zin)); % An array of roots is a column in Matlab conventions.

if(x(1) == 0)
  error('The highest order polynomial coefficient must not be 0.');
elseif(degree <= 0)
  error('Empty or constant polynomial.');
elseif(degree == 1)
  zout = zout - x(2)/x(1);
  return
elseif(degree == 2) % See Numerical Recipes section 5.6.
  rts = quadraticRoots(x);
  zout = zout + rts(2); % Initially set all output values to rts(2).
  % Go back and find if any input points were closer to rts(1).
  k = find(abs(zin-rts(1)) <= abs(zin-rts(2)));
  if(~isempty(k))
    zout(k) = rts(1);
  end
  return
end

if(nargin < 3) % Set the maximum number of iterations.
  maxNIter = 250;     
else
  maxNIter = fix(abs(maxNIter));
end

if(nargout > 1)
  zErr = zeros(size(zin));
  if(nargout > 2)
    nIter = zeros(size(zin));
  end
end

dx = (degree : -1 : 1) .* x(1 : degree);% dx=polyder(x) is much slower.

% BTW, there is a way to compute df in the same loop that computes f. However,
% for random coefficient polynomials, there is a difference in the 16'th digit between
% doing that and using polyder instead. For Wilkinson, there is a difference in the
% third digit. I do not have any evidence which is better, but my guess is that using
% polyder is more accurate.

for j=1:length(zin)
  if(DIAGNOSTICS & ~mod(j,100))
    fprintf(1,'\t\t\tPolishing root # %d\n', j);
  end
  z = zin(j);
  
  % Real polynomials have conjugate roots in the upper and lower
  % half plane. They pull equally on a position on the real axis.
  % Thus, doing a polish starting on the real axis takes a large
  % number of iterations. You want the offset from the real axis to
  % to be as large as possible without converging to a distant root.
  % Random coefficient polynomial roots always seem to be have a
  % separation of at least 1/(2*degree) so 1/(100*degree)*i seems
  % very safe.
  
  % The following is fine for random coefficient but not for 
  % x=poly([1 1e-10 1e-20]). If you start a polish at 1e-20 then 
  % z = z + 1/(100*degree)*i will move you quite far away and you
  % converge to 1e-10 instead of 1e-20.
  
  if(~imag(z)) % Move a real starting guess off the real axis.
     z = z + 1/(100*degree)*i;
  end

  % BTW, the Newton correction is f(z)/f'(z). Even if the denominator
  % is extremely small, do not arbitrarily set it to a larger value.
  % I have seen a naturally occurring polynomial where several of the
  % lowest order coefficients were on the order 1E-45. Thus the
  % polynomial and its derivative evaluated to 1E-45 at 0. The numerator 
  % and denominator were microscopic, but their ratio was not.
      
  NIter = 0;
  finished = 0; % Quit when finished == 2.
  while(finished < 2)
    NIter = NIter + 1;
    if(isinf(z) | isnan(z)) % Problems computing correction.
      temp = rand(1,1);
      z = cos(temp) + i*sin(temp); % Start at random position.
    else 
      f  = innerpart(x ,z);
      df = innerpart(dx,z);
      
      % Have we have found the root? This is not a perfect test.
      % If x=[1,1,zeros(1,14),-1,-1] then polyval(x,-1.00000000043683)==0 
      % but -1.00000000043683 is not a root. -1 is a double root.
      
      if(f == 0)
         correction = 0;
      elseif(~df)
         correction = Inf;
      elseif(abs(z) <= 1)
         correction = f/df;
      else
         correction = z*f/df;
      end
      
      z = z - correction;
      absZ = abs(z); % Used in the next line and a later section.
      relCorr = abs(correction)/absZ;
      
      % I did a test with one 1,000-degree and one 10,000-degree random coefficient
      % polynomial. Every root satisfied |final correction| < 1e-13*abs(z) and some
      % were .98e-13*abs(z). However, only 7% satisfied |final correction| < 1e-14*abs(z).
      % Therefore if |final correction| < 1e-13*abs(z) this is probably as good as it
      % gets but the following code ensures you do one more iteration. Perhaps you
      % can decrease from 9e-13 to 2e-13 with one more iteration. Note, you do not want 
      % a test like: if(absCorr<1e-16) because if the root equals 100, a correction
      % in the 16'th digit would be on the order of 1e-14. You need to bring abs(z)
      % into the picture.
      %
      % I have seen a 1,000 degree polynomial where the last 3 corrections for 1 root
      % were 1.1e-013 1.8e-016 and 6.9e-017. Surprisingly, you needed 2 corrections
      % after 1e-13. This is not converging quadratically. Also, amazingly, even when
      % the current correction was 1.8e-16, another iteration helped.
      
%       if(NIter >= maxNIter | relCorr < 1e-13)
      if(NIter >= maxNIter | abs(correction)<1e-14) % For degree 1e5 relCorr never got small enough.
         finished = finished + 1;
      end
      
      if(DIAGNOSTICS == 2)
         disp(['z=',num2str(z,8), ...
               '  r=',num2str(abs(z),8),'  theta=',num2str(angle(z)/pi*180,8), ...
              '   |correction| = ', num2str(abs(correction),6)]);
      end
    end
  end
      
  if(DIAGNOSTICS ~= 0)
     fprintf('\t\t\tNIter = %d\n\n',NIter);
  end
  
  if(nargout > 1)
     zErr(j) = relCorr;
     if(nargout > 2)
      nIter(j) = NIter;  
    end
  end
      
  % If real or imaginary part is close to 0, assume it is zero.
  if(absZ<=1)
     if(abs(imag(z)) < zErr(j))
        z = real(z);
     end
     if(abs(real(z)) < zErr(j))
        z = imag(z)*i;
     end
  else
     zInv = 1/z;
     if(abs(imag(zInv)) < zErr(j))
        z = real(z);
     end
     if(abs(real(zInv)) < zErr(j))
        z = imag(z)*i;
     end
  end
    
  % If, on the next line you use max(1E-7, 1/(5*degree)) and do the same
  % in the definition of gap in lroots() and and deflateFFT(), then you can
  % usually, but not always, factor Wilkinson. However, tlroots(1000,886) 
  % fails although the tlroots(1000, 1:885) works.
    
  % Even using abs(f)<1e-2 will not converge with Wilkinson because absCorr=.002.
  if(relCorr < 1e-3 & abs(f)<1e-13)
    zout(j) = z;
  else
    zout(j) = Inf;
 end
 
 zout(j) = z; % Debugging only
 
end

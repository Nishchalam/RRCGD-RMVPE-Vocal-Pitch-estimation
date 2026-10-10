% function [zout, zErr, nIter] = newton3(x, zin, maxNIter);
%
% PURPOSE:
%   Newton's method 3. If the current estimate, z, is outside the unit circle,
%   do z=1/z and polish against the flipped polynomial.
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
%   zErr - An error estimate. See comments below for an example where
%          the error estimate is too small by a factor of 1e38.
%
%          Also, if w is the Wilkinson polynomial, laguerre3(w,20.001) returns
%          an error estimate of 0 even though the actual error is 4e-5. However,
%          if you start at 20.01 or 20.0001 the error estimate is about 4e-5
%          which is correct. I.e., error estimate can be very wrong.
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
%
%
%   Actually we will not use polyval(). There are 2 natural ways to
%   evaluate a polynomial x at a value z.
%   1. f = polyval(x,z)
%   2. temp = filter(1, [1,-z], x)
%      f = temp(length(x)
%
%   But then 2. gives rise to a third and a fourth way.
%   3. Encapsulate 2. in a function that accepts a single value z.
%   4. Encapsulate 2. in a function that accepts an array of z values.
%
%   Timing tests using a 1023 degree random coefficient polynomial were
%   run using each of the 4 methods. Method 2. was the fastest. Method 3
%   was 5% slower. Method 4 was 10% slower. Method 1 was 15% slower.
%   Unfortunately, the fastest method is also the hardest to read.
%   Nonetheless I will use the fastest method because I need the speed.

% USES:
%   flip, normalizepow2, innerpart, quadraticRoots
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

% Also, [y,e]=newton3(w,15.001)
% y =   15.0061
% e =  2.0392e-040  % No this is not a good error estimate!

function [zout, zErr, nIter] = newton3(x, zin, maxNIter);

DIAGNOSTICS = 0; % 0=none, 1=minimal, 2=full diagnostic messages

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
  zErr = zeros(size(zout));
  return
elseif(degree == 2) % See Numerical Recipes section 5.6.
  rts = quadraticRoots(x); 
  zout = zout + rts(2); % Initially set all output values to rts(2).
  % Go back and find if any input points were closer to rts(1).
  k = find(abs(zin-rts(1)) <= abs(zin-rts(2)));
  if(~isempty(k))
    zout(k) = rts(1);
  end
  zErr = zeros(size(zout));
  return
elseif(degree == 3)
   % To be done later. cubic1.m & cubic2.m do not work well on poly(2 1e-10 1e-20)
   % or poly([1 2 1e-10]). Exponential gain might fix them.
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

dx = (degree : -1 : 1) .* x(1 : degree); % polyder(x) is much slower.
flippedX = flip(x);
dflippedX = (degree : -1 : 1) .* flippedX(1 : degree); % polyder(flipped_x)

% BTW, there is a way to compute df & d2f in the same loop that computes f. However,
% for random coefficient polynomials, there is a difference in the 16'th digit between
% doing that and using polyder instead. For Wilkinson, there is a difference in the
% third digit. I do not have any evidence which is better, but my guess is that using
% polyder is more accurate.

for j=1:length(zin)
  if(DIAGNOSTICS & ~mod(j,100))
    fprintf(1,'\t\t\tPolishing root # %d\n', j);
  end
  z = zin(j);
  
  stop=0;
  if(j == 218)
%     stop=1;
  end
  
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
  
  ef = 0; % Define ef so you can do abs(f)<=10*ef below.
  NIter = 0;
  finished = 0; % Quit when finished == 2.
  while(finished < 2)
    NIter = NIter + 1;
    if(isinf(z) | isnan(z)) % Problems computing correction.
      temp = rand(1, 1);
      z = cos(temp) + i*sin(temp); % Start at random position.
    else
      % It would be nice if we could say: df=innerpart(dx,z) independent of |z|>1 or not.
      % This fails because d(flip(x))~=flip(dx). Thus we need two cases below.
      
      % On another matter, if |z|>1 you would think innerpart(x,z)==innerpart(flippedX,1/z)
      % and therefore we could use innerpart(x,z) once outside the if-else. This is true for
      % innerpart.m but not innerpart.dll. Why? My C code does not accurately compute 1/z.
      if(abs(z)<=1)
         if(rem(NIter,5) & ~finished)
           f  = innerpart(x, z);
           df = innerpart(dx, z);
         else
           [f,ef]   = innerpart(x, z);
           [df,edf] = innerpart(dx, z);
         end
      else
         zz = 1/z;
         if(rem(NIter,5) & ~finished)
           f  = innerpart(flippedX, zz);
           df = innerpart(dflippedX, zz);
         else
           [f,ef]   = innerpart(flippedX, zz);
           [df,edf] = innerpart(dflippedX, zz);
         end
      end
      
      if(f == 0)
         correction = 0;
      elseif(df)
         correction = f/df;
      else
         correction = Inf;
      end
      
      if(stop)
         fprintf('%25.17g+%25.17gi\n',real(z),imag(z));
      end
      
      absZ = abs(z);
      if(absZ <= 1)
         z = z - correction;
         if(absZ)
            % Warning in Matlab 5.3,  abs(1e197 + 1e195i)==0.
            relCorr = abs(correction)/abs(z); % Relative correction.
         elseif(x(end)==0) % 0 is supposed to be a root.
            relCorr = 0;
         else
            relCorr = Inf;
         end
      else
         z = z/(1 - z*correction); % I.e. 1/(1/z-correction) which gives a slightly different answer.
         relCorr = abs(correction)/abs(zz); % Relative correction.
      end
      
      % I have seen a 1,000 degree polynomial where the last 3 absolute corrections for 1 root
      % were 1.1e-013 1.8e-016 and 6.9e-017. Surprisingly, you needed 2 corrections
      % after 1e-13. This is not converging quadratically. Also, amazingly, even when
      % the current correction was 1.8e-16, another iteration helped.
      
      if(NIter >= maxNIter | relCorr < 1e-13 | abs(f)<=10*ef)
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
  
  % Flipflag is no longer used.
  
  % How do we handle the error estimate if we are later going to do z = 1/z?
  % Suppose z is real and positive and the final correction, c, is accurate.
  % Then the correct answer lies between z-c and z+c. However, if flipFlag==1,
  % then our final answer will be 1/z and not z. Therefore, the range will
  % be 1/(z+c) to 1/abs(z-c). Our final answer of 1/z will be in error:
  % 1/z - 1/(z+c) = c/(z*(z+c)) or 1/(z-c) - 1/z = c/(z*(z-c)). Assuming
  % c was correct, the largest error in 1/z is c/(z*abs(z-c)).
  
  % This has only been tested on the Wilkinson polynomial. For the root 1,
  % c/(z*abs(z-c) and c were equally good. For 2 and 3, c was a better
  % estimate than c/(z*abs(z-c) which was too large. For 4 and greater,
  % c/(z*abs(z-c) was excellent, never being off by more than a factor of 10.
  % By contrast, c was too small by a factor of 100 to 4,000! The following
  % is a huge improvement on the error estimates for the Wilkinson polynomial.
  
  % There was a glaring exception. newton3(w,15.001) returns 15.006 however
  % the final correction near 1/15 was 1e-43 and the other was 1e-41. Why?
  % For the first 14 iterations, f(z) was about 1e-5 + a tiny imaginary value.
  % However, at iteration 15 and 16 f(z) was 0 + a tiny imaginary value and
  % the corrections were sub-atomic in size. Perhaps the polynomial evaluation
  % is the correct evaluation of a perturbed polynomial and 15.006 is a root
  % of that polynomial. Anyway, be aware that the error estimate in this case
  % was off by a factor of 1e-38! However, newton3(w,15.01) did not do this.
  
  % Also, the final correction is actually the correction to the next-to-last
  % root estimate, not the final root estimate. This might make a difference.
  
  f  = abs(f);
  df = abs(df);
  mx = max([f, df, ef, edf]);
  if(~isnan(mx) & ~isinf(mx) & mx)
     f = f/mx;   df = df/mx;  ef = ef/mx;  edf = edf/mx;
  end
    
  temp = df*abs(df-edf);
  if(~isnan(temp) & ~isinf(temp) & temp)
     zErr(j) = (f*edf+df*ef)/temp;
  else
     zErr(j) = relCorr; % I would rather use this than Inf. Hopefully df~=0.
  end
    
  if(nargout > 2)
    nIter(j) = NIter;  
  end
  
    % Real roots almost never polish to have imaginary part exactly 0.
  % If it is close to 0, assume it is zero.
  % BTW, I have seen a case where using 1E-10 failed.
  if(abs(imag(z)) < zErr(j))
    z = real(z);
  end
  
  if(abs(real(z)) < zErr(j)) % Warning. This turned 7.0711e-026 +7.0711e-026i into 7.0711e-026i
     z = imag(z)*i;
  end
  
  if(abs(f)<=10*ef)
    zout(j) = z;
  else
    zout(j) = Inf;
    fprintf('\n\t<<FYI - newton3() is returning Inf for root number j = %g>>\n\n',j);
    hwrite(['          -------- FYI newton3() is returning Inf for root number j = ',num2str(j)]);
  end
end

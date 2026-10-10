% function [zout, zErr, relCorr, nIter] = laguerre3(x, zin, maxNIter)
%
% PURPOSE:
%   Method 3 Laguerre
%
% REQUIRED INPUT:
%   x - polynomial coefficents. x(1) is the constant term.
%   zin - an array of starting points.
%
% OPTIONAL INPUT:
%   maxNIter - Default = 250. Maximum number of iterations.
%              I have seen simple examples where 200 is too small. Thus 250 might
%              also be too small. This parameter is important because of the
%              possibility of periodic orbits which never converge. They oscillate
%              from far outside the unit circle to far inside the unit circle.
%
% REQUIRED OUTPUT:
%   zout - A column vector of polished roots. Inf at locations where it failed to
%          polish in maxNIter steps.
%
% OPTIONAL OUTPUT:
%   zErr - Error estimates for the roots.
%
%          You would think you could use the final polishing correction and USUALLY
%          this is a good estimate. However, I have seen a final polishing correction
%          of 1e-41 for a root near the unit circle of a 1,000 degree random 
%          coefficient polynomial. Also, let w be the Wilkinson polynomial. In a 
%          previous version laguerre3(w,20.001) had a final polishing correction of
%          0 even though the actual error was about 4e-5.
%
%          Therefore a different error estimate is used. Once you have "converged"
%          the Newton and Laguerre corrections seem equivalent. The Newton correction
%          is -f(z)/f'(z). Function innerpart() has the option to return an estimate 
%          of the error in the polynomial evaluations. Since we can get error estimates
%          for f(z) and f'(z) we can compute an estimate of the error in -f(z)/f'(z).
%          This is used for the error estimate.
%   relCorr - The final polishing correction divided by the root magnitude.
%   nIter - The number of iterations actually used for each root.

% USES:
%   innerpart, normalizePow2, quadraticRoots
%
% WRITTEN BY:
%   James W. Fox, 10-15-2002 and patched many times after that.
% 
% CHANGES:
% 2/19/2006 If the computation of the Laguerre correction overflows, use the Newton correction instead.
% 2/23/2006 If it terminates from doing too many iterations, return Inf if the final relative polishing 
%           correction was not small. This may indicate you found a periodic cycle.

% This can fail in a bad way and I do not know what to do about it. Let
% x=[1 4 6 4 1].*[1 10 100 1000 10000]; % -10 is a root of multiplicity 4.
% Then laguerre2(x,-10.23) and laguerre3(x,10.13) return values with imaginary
% parts of .0008 & .001. These are larger than the error estimates so they are
% not set to 0. Factorization programs need to know if a root is real or not.

% Strange things can happen. Let w=poly(1:20) then polyval(flip(w),.05)==0.
% This is as it should be because 1/.05==20 is a root of w. However, 
% polyval(flip(w),.04999990026458411) == 8.970602038971265e-014
% polyval(flip(w),.04999990026458412) == 0
% polyval(flip(w),.04999990026458413) == 7.505107646466058e-014
% Thus laguerre3(flip(w),.04999990026458412) returns .04999990026458412 not .05.
% laguerre3(flip(w),.04999990026458413) is no better.

% I have seen an example where you start very near -1 and near several roots. 
% laguerre2 and laguerre3 converge to different nearby roots. laguerre3() converged to
% a root farther from the initial starting spot. abs(correction2)<abs(correction1)
% and polyval(x,z-correction2)<polyval(x,z-correction2) so correction2 is chosen.
% However using z-correction1 led to the closer root!

% This has a bug. Let x=[1 5 10 10 5 1]. laguerre3(x,-1.9) gives a terrible answer. 
% Internally z alternates between a good value and a terrible value. However, each time
% we have a good value, it is "corrected" and becomes a terrible value. Perhaps this is the 
% reason. The correction process involves H = GSquared - d2f/f; Note
% f = 2.7756e-016 -8.6027e-023i
% d2f = 8.8818e-016 -1.6200e-020i
% These values are likely highly sensitive and subject to error and the ratio amplifies error.

% DOCUMENTED CHANGES:
% 2/19/2006 If the computation of the Laguerre correction overflows, use the Newton correction instead.
% 2/23/2006 If it terminates from doing too many iterations, return Inf if the final relative polishing 
%           correction was not small. This may indicate you found a periodic cycle.

function [zout, zErr, relCorr, nIter] = laguerre3(x, zin, maxNIter)

global DIAGNOSTICS
if(isempty(DIAGNOSTICS))
  DIAGNOSTICS = 0;
end

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

x = normalizePow2(x(:).'); % Matlab convention is that a polynomial is a row vector.
degree = length(x) - 1;
zout = zeros(size(zin)); % An array of roots is a column in Matlab conventions.

if(degree == 0)
  error('A constant non-zero polynomial has no roots.')
elseif(degree == 1)
  zout = zout - x(2)/x(1);
  zErr = zeros(size(zout));
  relCorr = zErr;
  return
elseif(degree == 2) % See Numerical Recipes section 5.6.
  rts = quadraticRoots(x);
  zout = zout + rts(2); % Initially set all output values to rts(2).
  % Go back and find if any input points were closer to rts(1).
  k = find(abs(zin-rts(1)) <= abs(zin-rts(2)));
  if(~isempty(k))
    zout(k) = rts(1);
  end
  zErr = zeros(size(zout)); % Just because you have a formula does not mean this is true.
  relCorr = zErr;
  return
end

if(nargin < 3)
  maxNIter = 250;  % Maximum number of iterations.
else
  maxNIter = fix(abs(maxNIter));
end

if(nargout > 1)
  zErr = zeros(size(zin));
  if(nargout > 2)
     relCorr = zeros(size(zin));
     if(nargout > 3)
        nIter = zeros(size(zin));
     end
  end
end

% We have already taken care of degree==1 or 2 above so degree >= 3.
dx = (degree : -1 : 1) .* x(1 : degree); % polyder(x) is much slower.
d2x = (degree-1 : -1 : 1) .* dx(1 : degree-1); % polyder(polyder(x))
flippedX   = flip(x);
dflippedX  = (degree : -1 : 1) .* flippedX(1 : degree); % polyder(flippedX)
d2flippedX = (degree-1 : -1 : 1) .* dflippedX(1 : degree-1);
% BTW, there is a way to compute df & d2f in the same loop that computes f. However,
% for random coefficient polynomials, there is a difference in the 16'th digit between
% doing that and using polyder instead. For Wilkinson, there is a difference in the
% third digit. I do not have any other evidence which is better, but my guess is
% that using polyder is more accurate.

for j=1:length(zin)
  if(DIAGNOSTICS & ~mod(j,1000))
    % BATCH MODE: disabled diagnostic output
    % fprintf(1,'\t\t\tPolishing root # %d\n', j);
  end
  
  % BTW, do not attempt to reset tiny values. I have seen
  % a naturally occurring polynomial where several of the lowest
  % order coefficients were on the order 1E-45. Thus the polynomial
  % and its derivative evaluated to 1E-45 at 0. Thus a polish that
  % started at 0 required f(z) and f'(z) to be about 1E-45.

  z = zin(j);
  ef = 0; % Define ef so you can do abs(f)<=10*ef below.
  NIter = 0;
  finished = 0; % Quit when finished == 2.
  while(finished < 2)
    NIter = NIter + 1;
    if(isinf(z) | isnan(z)) % Problems computing correction.
      temp = pi*rand(1, 1);
      z = cos(temp) + i*sin(temp); % Start at random position.
    else
      if(abs(z) <= 1)
        if(rem(NIter,5) & ~finished)
           f   = innerpart(x, z);
           df  = innerpart(dx, z);
           d2f = innerpart(d2x, z);
        else
           [f,ef]     = innerpart(x, z);
           [df,edf]   = innerpart(dx, z);
           [d2f,ed2f] = innerpart(d2x, z);
        end
      else
         zz = 1/z;
         if(rem(NIter,5) & ~finished)
            f   = innerpart(flippedX  , zz);
            df  = innerpart(dflippedX , zz);
            d2f = innerpart(d2flippedX, zz);
         else
           [f,ef]     = innerpart(flippedX, zz);
           [df,edf]   = innerpart(dflippedX, zz);
           [d2f,ed2f] = innerpart(d2flippedX, zz);
         end
      end
      
      % Have we have found the root? This is not a perfect test.
      % If x=[1,1,zeros(1,14),-1,-1] then polyval(x,-1.00000000043683)==0 
      % but -1.00000000043683 is not a root. -1 is a double root.
      if(f == 0) % We have supposedly found the root.
         correction = 0; % I have seen f==0 and ef==0 when z was not EXACTLY a root.
      else
        % x = [1,1,zeros(1,114),-1,-1]; % This has -1 as a double root.
        % z = -1.000000000812914 - 7.0334588403050095e-157*i; % This is virtually a root.
        % I did laguerre3(x,z,200) and on the first iteration:
        % f == 0 +1.326482327954642e-163i % f is virtually zero at z.
        % df == 1.885960266845288e-007 +1.631762219479584e-154i
        % G == 1.230142449010735e+009 -1.421775644575173e+156i % HUGE imaginary part.
        % G*G == -Inf -3.497973146723039e+165i
        % H == -Inf -1.748986653816523e+165i
        % temp == NaN - NaNi
        % The newly added section: if(isnan(temp) & abs(f)<1e-10) ... takes care of this,
        % at least for this polynomial.
        
        G = df/f; % This looks dangerous if f==1e-16.
        GSquared = G*G; % This can overflow or underflow.
        H = GSquared - d2f/f;
        temp = sqrt((degree-1)*(degree*H - GSquared));
        
        if(isnan(temp)) % Laguerre correction overflowed. Use Newton correction instead. 2/19/2006
           correction1 = f/df;
           correction2 = correction1;
        else
           if(abs(G) ~= abs(temp))
             correction1 = degree/(G + temp);
             correction2 = degree/(G - temp);
           else
              if(z==0)
                 correction1 = -1e-25*i; % This is completely ad hoc. It works for z=0; &
                 correction2 = correction1; % x=conv([-1e-90 0 0 1] , [1 0 0 0 1e-100])
              else
                 correction1 = Inf;
                 correction2 = Inf;
              end
           end
        end
        
        if(abs(correction1) <= abs(correction2))
          correction = correction1;
        else
          correction = correction2;
        end
      end
      
      absZ = abs(z); % Used in the next section and a later section.
      if(absZ <= 1)
         z = z - correction;
         if(absZ)
            % Warning in Matlab 5.3,  abs(1e197 + 1e195i)==0.
            relCorr(j) = abs(correction)/absZ; % Relative correction.
         elseif(~x(end)) % 0 is supposed to be a root.
            relCorr(j) = 0;
         else
            relCorr(j) = Inf;
         end
      else
         temp = z*correction;
         if(1 ~= temp)
            z = z/(1 - temp); % Equivalent to 1/(1/z-correction).
            % Let r=the exact root and let c=correction computed above.
            % Then c is the correction for 1/z. I.e., |1/z-1/r|==c so
            % |z-r|==c*|rz|. We do not know r but we assume it is close to
            % z so |z-r|, the error in z, is approximately c*|z|^2. Thus
            % the relative error is c*|z|.
            relCorr(j) = abs(correction)*absZ;
         else
            z = Inf;
            df = 0;
            ef = 0;
            f = Inf;
         end
      end
      
      % I have seen a 1,000 degree polynomial where the last 3 corrections for 1 root
      % were 1.1e-013 1.8e-016 and 6.9e-017. Surprisingly, you needed 2 corrections
      % after 1e-13. This is not converging quadratically. Also, amazingly, even when
      % the current correction was 1.8e-16, another iteration helped.
      
      if(NIter >= maxNIter | relCorr(j) < 1e-15 | abs(f)<=10*ef)
         finished = finished + 1;
      end
      
      if(DIAGNOSTICS == 2)
         % BATCH MODE: disabled diagnostic output
         % disp(['z=',num2str(z,8), ...
         %       '  r=',num2str(abs(z),8),'  theta=',num2str(angle(z)/pi*180,8), ...
         %      '   |correction| = ', num2str(abs(correction),6)]);
      end
    end
  end
  
  % y=sin(0:2*pi/100:2*pi);
  % a=laguerre3(y,0)            % a =-3.9007e-015
  % d=deconv(y,[1,-a]);         % Deflate that root
  % b = [.6422739514522119 -.4929268048301818i, .9798339281043483 +.7519943886134968i]
  % Then b is a periodic cycle of length 2 for d. You finish because NIter=maxNIter+1 but you never found a root.
  if(NIter >= maxNIter & relCorr(j) > 1e-7) % Probably a periodic cycle  2/23/2006
     z = Inf;
  end
     
  if(DIAGNOSTICS == 2)
     % BATCH MODE: disabled diagnostic output
     % fprintf('\t\t\tThe number of iterations  = %d\n\n',NIter);
  end
  
  % ef and edf are error estimates for f and df. Thus we can obtain an estimate
  % for the Newton correction -f/df and use this for an error estimate.
  % 
  % Physics courses teach about error propagation in computed measurements.
  % For example suppose you measure two quantities x and y and you estimate
  % the errors in the measurements as dx and dy. If you then compute x+y, or x/y
  % what is the error in these computations? Numerous sites on the web will tell
  % you that the error in x+y is sqrt(dx^2+dy^2) and the error in x/y is 
  % (x/y)*sqrt((dx/x)^2 + (dy/y)^2). This is true provided you make an assumption
  % which is true in physics labs - you have taken several measurements of x and y.
  % You computed averages and use them for x and y and you computed standard 
  % deviations of the errors and used them for dx and dy. There is, however, an 
  % alternative - use the mean, not standard deviation, of the absolute values of 
  % the deviations for dx and dy. Statistically this is a poorer procedure if you 
  % have taken several measurements. However, what if you have only taken one 
  % measurement, which is our case in using ef and edf as error estimates. In that
  % case using an alternative formula may provide a larger and safer estimate for 
  % the error in the computed value. A good site that discusses the two methods is
  % http://www.rit.edu/~vwlsps/uncertainties/Uncertainties.html
  
  % For example consider z=x+y. Suppose we have errors dx and dy. I.e. if we had
  % done the measurement correctly we would have gotten x+dx and y+dy. Then we
  % would have gotten (x+dx)+(y+dy) instead of x+y. The error is thus dx+dy. Either
  % measurement could be negative (or complex) so a safe error estimate is
  % |dx| + |dy|. This is larger than the standard physics estimate of 
  % sqrt(dx^2+dy^2). To see this, square both; the first is larger by 2*|dx|*|dy|.
  % Since dx and dy are small, the standard physics error estimate says their product
  % is so small it can be ignored. However, if you have only taken one measurement
  % and you want a safer error estimate, then you should use |dx| + |dy|.
  
  % How do we estimate the error in x/y assuming possible errors dx and dy. The correct
  % value is (x+dx)/(y+dy). If you use long polynomial division, the first equals 
  % x/y + (y*dx - x*dy)/(y*(y+dy)). Therefore the error in using x/y is 
  % (y*dx - x*dy)/(y*(y+dy)). You get the same value if you use algebra to simplify
  % (x+dx)/(y+dy) - x/y. Some of these values could be negative or complex so a safer
  % error estimate is |(y*dx - x*dy)/(y*(y+dy))| <= (|y*dx|+|x*dy|)/(|y| * | |y|-|dy| |).
  % Note, |y+dy| >= |y|-|dy| is the backwards triangle inequality.
  
  % Here is a numerical example. Suppose x=16 dx=0.2 y=8 dy=-0.1. The exact error is
  % (x+dx)/(y+dy)-x/y = 0.05063291139240. Using my formula gives 0.05063291139241
  % which is off by 1 in the 14th decimal place. The standard formula from physics
  % using sqrt gives 0.035 which needs to be increased by 66%! My formula is safer.
  
  % Let x=conv([1 0 0 0 0 0 0 0 0 0 0 0 -1e-300],[1 0 0 0 0 0 0 0 0 0 0 0 -1e300]);
  % In a previous version laguerre3(x,0)==0 because temp, below, underflowed. f, df etc
  % are TINY but not 0, however their product was zero. If you rescale f, df, ef, & edf
  % and then compute the ratio, everything is ok.
  
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
     zErr(j) = relCorr(j); % I would rather use this than Inf. Hopefully df~=0.
  end
    
  if(nargout > 2)
    nIter(j) = NIter;  
  end
  
  % Recall the line above, relCorr(j) = abs(correction)*absZ; and the discussion
  % before it. The same argument applies to zErr. I.e., if |z|>1 then we need
  % to multiply the computed error for 1/z by |z|^2.
  %
  % For some reason, in the Wilkinson polynomial, using |z|^2 gives error estimates
  % that are too large and it would be much better to use |z|. However I think the
  % theory is correct from the following experiment which was done several times.
  % z=1e6+rand(1,1); t=rand(1,1000); t=conv(t,[1,-z]); [a,b]=laguerre3(t,z+100); 
  % Then a-z is the exact error. It was always either 0 or about 2e-10. b==6e-10 which
  % is a great error estimate. Note, z has 6 digits to the left of the decimal and
  % 10 digits to the right of the decimal. Thus the error must be at least 1e-10. 
  % If you do not do the following section then zErr==6e-16 which can not be correct.
  
  % zErr is no longer used. I switched to using relCorr.
  
  if(absZ > 1)
     zErr(j) = zErr(j)*absZ*absZ;
  end
  zErr(j) = max(zErr(j), 5*eps*absZ); % Totally ad hoc.
  
  % Zero out the real or imaginary part if it is extremely small.
  if(abs(imag(z)) < 2*absZ*min(1e-12,relCorr(j)))
     z = real(z);
  end
  
  % z=1e10+1e30i; x=poly([1e10 z conj(z)]); laguerre3(x,1e30i)== 0+1e30i. The real part 
  % was 9.999999999999996e9 until the following lines makes it 0 which is wrong. However,
  % in another poly it takes 1.52e-065 - 1.e-020i and correctly turns it into -1e-20.
  % Obviously, I do not really know what to do.
  if(abs(real(z)) < 2*absZ*min(1e-12,relCorr(j)))
    z = imag(z)*i;
  end
  
% I do not know what to do. Let x=poly([1e6 1e40 1e80]); Then z=9.999999999999999e+005,
% f=1e-22, & ef=2e-37. I do not know if it is significant but ef/f=2e-15. It is
% not true that abs(f)<10*ef so the original code set z to Inf even though z is 
% almost perfectly accurate. For the time being I am going to comment out these lines
% and just accept z.

  zout(j) = z;
  
%  if(abs(f)<=10*ef)
%    zout(j) = z;
%  else
%    zout(j) = Inf;
%    fprintf('\n\t<<FYI - laguerre3() is returning Inf for root number j = %g>>\n\n',j);
%    hwrite(['<<          -------- FYI laguerre3() is returning Inf for root number j = >> ',num2str(j)]);
%  end
end

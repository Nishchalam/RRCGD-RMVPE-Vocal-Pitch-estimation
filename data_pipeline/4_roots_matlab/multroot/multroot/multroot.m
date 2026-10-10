% function [z1d, z2d, parms] = multroot(p, printFlag, tol, thresh, growf)
%
%  MultRoot calculates roots and multiplicities of a polynomial p,
%        p(x) = a_1 x^n + a_2 x^n-1 + ... + a_n x + a_n+1,
%  given as an (n+1)-row vector p = (a_1, ..., a_n+1) 
%
%  MultRoot was developed by Zhonggang Zeng  (email: zzeng@neiu.edu)
%  This code is freely released for research exchange only. The author 
%  is not responsible for any damage from using this code.
%
%  INPUT :       p -- polynomial coefficients.
%        printFlag -- Default=1. Set to 0 to suppress printing of results (optional).
%              tol -- Default=1e-10. Initial residual tolerance (optional).
%           thresh -- Default=100*tol. Zero singular value threshold (optional).
%            growf -- Default=100. Residual growth factor (optional).
%
%  OUTPUT :
%            z1d -- a 1d array of polynomial roots where a root with
%                   multiplicity k occurs k times.
%            z2d -- roots of polynomial p and multiplicities:
%                      z2d(:,1) = distinct roots
%                      z2d(:,2) = corresponding multiplicities
%            parms -- [forward error; backward error; pejorative condition number]
%            
%  CALL :   The simplest way to call:
%             >> z1d = multroot(p);
%           For expert users, 1-4 parameters can be entered
%           For example:
%             >> [z1d,z2d,parms] = multroot(p, 0, 1e-12,1e-9,10);
%           etc. 

% JWF
% Modifications:
%  Zhonggang Zeng, November 2003, ????
%
%  James W. Fox, April-May 2004:
%    Corrected several minor bugs - some low degree polynomials could cause an abort.
%    Added 3 new arguments - printFlag, z1d, and parms. Removed 3 arguments.
%    Added new functions z2d_z1d, expGain, displayRoots.
% JWF   
%function [z,f_err, b_err, cond] = ...
%                   multroot( p, tol, thresh, growf )
function [z1d, z2d, parms] = multroot( p, printFlag, tol, thresh, growf )

% JWF
% multroot(['abc']) produces an answer when it should not.
% multroot([1 Inf 1]) and multroot([1 NaN 1]) have ugly aborts.
% The new code fixes these problems.
if(ischar(p))
      error('$1 must be a real array.');
   elseif(max(isnan(p)))
      error('$1 must not contain NaN or -NaN.');
   elseif(max(isinf(p)))
      error('$1 must not contain Inf or -Inf.');
   end
   
% JWF
% Added new parameter, printFlag, and chaged gamma to growf
%
% set default parameters if not provided
%   if nargin == 1
%       tol = 1.e-10; thresh = tol*100; gamma = 100;
%   elseif nargin == 2
%       thresh = tol*100; gamma = 100;
%   elseif nargin == 3
%       gamma = 100;
%   else
%       gamma = growf;
%   end;
   
% Set default parameters if not provided
   if nargin == 1
       printFlag = 1; tol = 1.e-10; thresh = tol*100; growf = 100;
   elseif nargin == 2
       tol = 1.e-10; thresh = tol*100; growf = 100;
   elseif nargin == 3
       thresh = tol*100; growf = 100;
   elseif nargin == 4
       growf = 100;
   end;
   
% JWF
% Useful error checks. You could do more on the reasonableness of tol, thresh, & growf.
   if(ischar(printFlag) | length(printFlag)~=1)
      error('$2 must be a real constant.')
   elseif(ischar(tol) | length(tol)~=1 | tol <= 0)
      error('$3 should be a small positive real number, eg. 1e-10.');
   elseif(ischar(thresh) | length(thresh)~=1 | thresh <= 0)
      error('$4 should be a small positive real number, eg. 1e-8.');
   elseif(ischar(growf) | length(growf)~=1 | growf <= 0)
      error('$5 should be a relatively large positive real number, eg. 100.');
   end;
   
   %
   % transpose if p is a column vector
   %
   
% JWF
% Check that p is not mxn where m and n > 1. multroot(rand(3,3)) has an ugly abort.
%  [m,n] = size(p);            
%  if m > n, p = p.'; end;
   [m,n] = size(p);
   if(m ~= 1 & n ~= 1)
      error('$1 must either be a row or a column vector.');
   elseif(n == 1)
      p = p.'; % Convert column vector to row vector.
   end

   %
   % clear leading/trailing zeros
   %
   n = length(p);
   q = p;
   if q(1) == 0 | q(n) == 0 
       jj = find(p); 
       j1 = min(jj); j2 = max(jj);
       q = p(j1:j2);
% JWF 
% multroot([0 0]) has an ugly abort because q == [] so q/q(1) aborts.
%      q = q/q(1);
       if(~isempty(q))
          q = q/q(1);
       end
   else
       j1 = 1; j2 = n;
   end;
   
% JWF
% multroot([1 0 0] has an ugly abort because q==1 which has no roots.
% Also print a message if it was the zero polynomial.   
   degree = length(p)-j1;
   if(isempty(q)) % p == [0 0 ... 0 0]
      error('Every complex number is a root of the zero polynomial.');
   elseif(length(q) == 1) % p = [1 0 0 0 ... 0] etc.
      if(degree)
         z1d = zeros(degree,1);
         z2d = [0, degree];
         parms = [0 0 0];
         return;
      else
         error('A non-zero zero-degree polynomial does not have any roots.');
      end
   end
      
   %
   % exponential scaling
   %
% JWF
   if(length(q) >= 4 & max(abs(q(end-2:end))) < 1e-10*max(abs(q(1:end-3))))
      % Initially I commented out the following scaling and instead used expGain.
      % Zeng's exponential scaling only used the first and last coefficient.
      % However, when expGain was applied to x=[1 1e-200 0 0 -1e-80] the result was
      % a single root of multiplicity 4 instead of 1e-20*(the fourth roots of unity)
      % while Zeng's procedure would get the desired result. Why? expGain tried to fit
      % all the coefficients including 1e-200. It came up with a gain factor of 1e-7.
      % When this was applied, the resulting polynomial had roots approximately equal to
      % 1e-13*(the fourth roots of unity). Zeng's gcdroot decided this was a single root
      % of multiplicity 4. By contrast, if you use Zeng's original scaling using only
      % the first and last coefficient, the gain factor is 1e20. The resulting polynomial
      % is the fourth roots of unity. Gcdroot could easily handle that. 
      %
      % In other words, gcdroot will have problems if there are several microscopic roots.
      % If there are, you need to apply a scaling which greatly increases their magnitude.
      % If there are 2 or more really tiny roots then the last 3 coefficients will be really
      % tiny. The above test uses this and if it decides there are probably really tiny
      % roots, then it uses Zeng's original scaling. If not, it uses expGain. This at least
      % works for the above polynomial.
      %
      % This is Zeng's original multroot scaling
      m = length(q)-1;
      %k = ceil(log(abs(q(m+1)))/(m*log(2)));
      k =       log(abs(q(m+1)))/(m*log(2));
      c = 2^(-k);  %c = 1;
      q = q.*(c.^[0:m]);
      q = q/q(1); % Some things, like bkerr in gcdroot, need q(1) to be 1
      gain=1/c;
   else
      [q,gain] = expGain(q);
      q = q/q(1); % Some things, like bkerr in gcdroot, seem to need q to be monic.
   end
   
% JWF
% This is only because I renamed gamma to growf.   
%  [y,bke] = gcdroot(q, tol, thresh, gamma);
   [y,bke] = gcdroot(q, tol, thresh, growf);
   
   if bke < 1.0d-2 & max(y(:,2)) > 1
% JWF
% This is only because I renamed z to z2d.      
%      [z,bkerr,pjcnd,job] = pejroot(q,y(:,1).',y(:,2)');
       [z2d,bkerr,pjcnd,job] = pejroot(q,y(:,1).',y(:,2)');
       
       if job == 1
           if j2 < n
% JWF
% This is only because I renamed z to z2d.      
%             z = [z;0,n-j2];
              z2d = [z2d;0,n-j2];
           end;
           
% JWF
% This is only because I changed c to 1/gain. See expGain above.
%          z(:,1) = z(:,1)/c;
           z2d(:,1) = gain*z2d(:,1); % Remove the effects of exponential scaling
           pjcnd = spcond(z2d);
       else
% JWF
% p does not have exponential scaling. q does so it is better to work with q. Plus there is
% a second problem if you use p instead of q. Suppose x=[0 1 1].  The first two lines of
% rtsinfo are
%
% n = length(p)-1;
% z = [roots(p),ones(n,1)];
%
% In this case n==2 but roots(p) has only one root. Thus you abort when you try to
% concatenate roots(p) and ones(n,1). 
%          [z2d,bkerr,pjcnd] = rtsinfo(p);
           [z2d,bkerr,pjcnd] = rtsinfo(q);
% JWF       
           if j2 < n
              z2d = [z2d;0,n-j2];
           end           
           z2d(:,1) = gain*z2d(:,1); % Remove the effects of exponential scaling
       end;
   else
% JWF
% See the comments above about why it is better to use q and not p.
%      [z2d,bkerr,pjcnd] = rtsinfo(p);
       [z2d,bkerr,pjcnd] = rtsinfo(q);
% JWF       
       if j2 < n
          z2d = [z2d;0,n-j2];
       end           
       z2d(:,1) = gain*z2d(:,1); % Remove the effects of exponential scaling
   end;
   ferr = 2*bkerr*pjcnd;
   
   z2d = displayRoots(z2d,0); % Clip small parts and sort; do not print yet.
   z1d = z2d_z1d(z2d); % Create 1d array from 2d array.
   
   if(nargout > 2)
      parms = [ferr; bkerr; pjcnd];
   end
   
   %
   % show off results
   %
   
% JWF. 
% There is a new input variable, printFlag.
   if(printFlag)
      % BATCH MODE: disabled diagnostic output
      % fprintf('\n');
      % fprintf('THE CONDITION NUMBER:                         %g \n',pjcnd);
      % fprintf('THE BACKWARD ERROR:                    %6.2e \n',bkerr);
      % fprintf('THE ESTIMATED FORWARD ROOT ERROR:      %6.2e \n',...
      %                                                 ferr);
      % fprintf('\n');
      % displayRoots(z2d,1);
   end
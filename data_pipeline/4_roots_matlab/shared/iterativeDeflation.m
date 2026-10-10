% function [z, zErr, rc, outParms] = iterativeDeflation(x, z, zErr, ID_NUM)
%
% PURPOSE
%   Given an incomplete but polished set of upper half-plane roots
%   of a real polynomial, try to find the rest of the upper half-plane
%   roots. It is ok to have a double root in z, but it had better be 
%   a genuine double root; i.e., you probably want to use uniq()
%   immediately before calling this routine.
%
%   Deflate roots that are not too close to the unit circle in the
%   Fourier domain. Deflate unit circle-ish roots in the time domain.
%
% REQUIRED INPUT
%   x - The polynomial.
%   z - An incomplete set of roots in the upper half plane.
%       They MUST already be polished and they MUST be unique.
%   zErr - An array or error estimates for z. If [] is supplied, an array of 1e-7's
%          is used. If only a single number is supplied, it is repeated. 
%   ID_NUM - A number that identifies this polynomial. It is used in
%            information messages to the history file. If this is 0,
%            messages are suppressed.
%
% REQUIRED OUPUT
%   z - A larger, hopefully complete, set of roots. However it may be
%       imcomplete; check the return code, rc.
%
% OPTIONAL OUPUT
%   rc - The return code.
%        rc < 0 indicates failure. 
%           -1 Too many roots. Probable polish error.
%           -2 Correct degree, but unfactor() did not confirm.
%           -4 Deflation provided no new roots.
%           -5 Correct degree, but unfactor() did not confirm and no new roots.
%           -6 Degree too small and no new roots.
%        rc == 0 means polynomial division was not needed; all roots were input.
%        rc>0 is the number of polynomial divisions before all roots were found.
%   outParms - an array of output parameters.
%              outParms(1) = the number of roots obtained by deflation.
%              outParms(2) = an estimation of max(abs(|x - unfactor(z)|)). However,
%                            if it was necessary to use deflation to find roots,
%                            it may be several orders of magnitude too large.

% USES
%   global DEFLFILE.
%   hpdegree, hwrite, ixFirstNZ, ixLastNZ, nonPolish, normalizePow2,
%   polycmp, laguerre3, unfactor, uniq, updateR.
%
% AUTHOR
%   James W. Fox 4/22/2001 and patched several times after that.

% By the way, this routine uses unfactor(z,1) to unfactor a set of roots. If you are
% only going to factor low degree you would be better off using poly(leja(z)) in this code.
% Poly(leja(z)) will probably be more accurate, but significantly slower for large degree.
% I tried leja(z) for a million degree polynomial and I killed it after 60 hours of 
% uninterrupted computation. By contrast, unfactor(z,1) did a van der Corput sort and 
% finishes the unfactorization in 32 hours. 

function [z, zErr, rc, outParms] = iterativeDeflation(x, z, zErr, ID_NUM)
global DEFLFILE; % Set in RootsIni

x=x(:)';
if(any(imag(z) < 0))
   error('All roots, $2, must be in the upper half plane.');
end
z = uhp(z); % It will error off if there are lower half plane roots.
z = z(:);
zErr = zErr(:);

if(isempty(x))
   z = [];
   zErr = [];
   rc = 0;
   outParms = [0 0];
   return
elseif(prod(size(z)) ~= prod(size(zErr)))
   error('Must have size($2)==size($3)');
elseif(isempty(zErr))
   zErr = 1e-7 + zeros(size(z));
elseif(length(zErr == 1))
   zErr = abs(zErr) + zeros(size(z));
else                          
  if(length(zErr) ~= length(z))
    error('length($3) must equal length($4).')
  elseif (min(zErr) < 0)
    error('$4 must be not contain negative values.')
  end
end

%---- Remove leading and trailing zeros and normalize ----
z = z(find(z ~= 0)); % Remove any zeros in z
ixLNZ = ixLastNZ(x);
NZeros = length(x) - ixLNZ; % Multiplicity of 0.0 as a root. Used at the very end.
x = normalizePow2(x(ixFirstNZ(x) : ixLNZ));
Lx = length(x);
degree = Lx - 1; % Note, degree == 0 is possible. Consider x=[0 0 1 0]
Lsp = 2^nextpow2(Lx);
sp_x = fft(x,Lsp); % The spectrum of x.
rc = 0;   % Initialize the return code
nDeflations = 0;
maxErr = Inf;
max_error_allowed = 3e-4; % Success if: |x-unfactor(z)| < max_error_allowed
% max_error_allowed should depend on zErr but currently zErr is not very accurate.
deflatedR = []; % Root radii obtained by deflation. I write them to a file.
deflatedRErr = []; % Root radii obtained by deflation. I write them to a file.

% We want to do deflation in the Fourier domain. However, we can not deflate a
% root in the Fourier domain if the spectrum of that root is too small; i.e. if
% it is near an n'th root of unity. Therefore, split the roots into roots that
% are extremely close to the unit circle, zUC, and those that are not, z. 
% My test for detecting unit circle roots is a wild guess, based only on the
% fact that for 150,000 degree random polynomials, you can get one root this
% close to the unit circle by chance. Also this code worked to factor 150,000.

% I wish the following test was more intelligent.
test = abs(abs(z)-1) < 1e-12; % True/false where distance to unit circle is small/large.
ixUC    = find(test);  % Indicies of roots close to the unit circle.
ixFar   = find(~test); % Indicies of roots far from the unit circle.
zUC     = z(ixUC);
zFar    = z(ixFar);
zUCErr  = zErr(ixUC);
zFarErr = zErr(ixFar);
clear test; clear z; clear zErr;

if(isempty(zUC)) % Compute the spectrum of roots close to the unit circle.
  sp_zUC = 1 + zeros(1, Lsp);
else   
  [junk, sp_zUC]= unfactor(zUC, 1, 0, Lsp);
end

if(isempty(zFar)) % Compute the spectrum of roots far from the unit circle.
  sp_zFar = 1 + zeros(1, Lsp);
else  
  [junk, sp_zFar] = unfactor(zFar, 1, 0, Lsp);
end

done = 0;
while(~done)
  hpdegzFar = hpdegree(zFar); 
  hpdegUC   = hpdegree(zUC); 
  hpdeg     = hpdegzFar + hpdegUC;
  
  if(hpdeg > degree) % Give up.
    z = [zFar ; zUC];
    zErr = [zFarErr ; zUCErr];
    done = 1;
    rc = -1;
  elseif(hpdeg == degree)
    u = ifft(sp_zFar .* sp_zUC);
    u = real(u(1:Lx));
    z = [zFar ; zUC];
    zErr = [zFarErr ; zUCErr];
    done = 1;
    maxErr = polycmp(x,u);
    if(maxErr < max_error_allowed)
       rc = nDeflations;
    else
       rc = -2;
    end
  else % hpdeg < degree
    nDeflations = nDeflations + 1; % The number of deflation steps.
    if(hpdegzFar == 0)
      quo = x;
    else
      if(any(sp_zFar == 0)) % I have seen sp_zFar underflow to 0 at 1 frequency.
         t_zFar = real(ifft(sp_zFar, Lsp));
         t_zFar = t_zFar(1: (hpdegzFar+1));
         quo = deconv(x, t_zFar); % Deflation in the time domain.
      else
         quo = ifft(sp_x ./ sp_zFar); % Deflation in the Fourier domain.
         quo = real(quo(1: degree - hpdegzFar + 1)); 
      end
    end
    
    if(~isempty(zUC))
      xUC = ifft(sp_zUC);
      xUC = real(xUC(1:(hpdegUC + 1)));
      quo = deconv(quo, xUC);
    end
    quo = normalizePow2(quo); 
    c = roots(quo);
    c = c(find(imag(c) >= 0)); % Upper half-plane roots.
    c = laguerre3(quo, c);     % Polish against quo.
    c = nonPolish(c, []);      % Polish against quo.
    [c,cErr] = laguerre3(x, c);    % Polish against x.
    [c,cErr] = nonPolish(c, cErr); % Remove roots that polished to 0, Inf, or lower half plane.
    
    if(isempty(c))  % nonPolish() removed everything.
      z = [zFar;zUC];
      zErr = [zFarErr; zUCErr];
      done = 1;
      rc = -4;
    else
      % Keep track of root radii obtained by deflation -deflateR. Write 
      % them to a file for later study. Doing a histogram of these
      % roots shows where you need to make the grid smaller. I.e.,
      % If the histogram shows a huge number of roots obtained
      % by deflation with radii near 1-(5/degree), then you need to
      % make your grid smaller there. Note, c may contain a false
      % double root, but that is irrelevant for this purpose.
      
      temp = abs(c); % Radii gotten by deflation
      j = find(temp > 1);  % Convert >1 radii to equivalent <1.
      if(~isempty(j))
        temp(j) = 1 ./ temp(j);
      end
      deflatedR = [deflatedR; temp];
      deflatedRErr = [deflatedRErr; cErr];
      [deflatedR,deflatedRErr] = uniq(deflatedR, 10*deflatedRErr);
      deflatedRErr = deflatedRErr/10; % Must be the same factor as in $2 on previous line.
      
      % Once, c obtained from the intermediate quotient had 2 real roots.
      % When they were polished against the original polynomial, one
      % polished to a previously discovered root and the other polished
      % to a complex root. This meant hpdeg>degree below.
      hpdeg = hpdeg + hpdegree(c);
      if(hpdeg > degree)
        [zFar,zFarErr,zUC,zUCErr,sp_zFar,sp_zUC,LNew]=...
                         updateR(c,cErr,zFar,zFarErr,zUC,zUCErr,sp_zFar,sp_zUC);
        if(~LNew)
          z = [zFar;zUC];
          zErr = [zFarErr; zUCErr];
          done = 1;
          if(hpdegree(z) > degree)
             rc = -1;
          else
             rc = -4;
          end
        end
      elseif(hpdeg == degree)
        [junk, sp_C] = unfactor(c, 1, 0, Lsp);
        u = ifft(sp_zFar .* sp_zUC .* sp_C);
        u = real(u(1:Lx));
        maxErr = polycmp(x,u);
        if(maxErr < max_error_allowed)
          z = [zFar; zUC; c];
          zErr = [zFarErr; zUCErr; cErr];
          done = 1;
          rc = nDeflations;
        else
          [zFar,zFarErr,zUC,zUCErr,sp_zFar,sp_zUC,LNew]=...
                         updateR(c,cErr,zFar,zFarErr,zUC,zUCErr,sp_zFar,sp_zUC);
          if(~LNew) % c contained no new roots.
            z = [zFar ; zUC];
            zErr = [zFarErr; zUCErr];
            done = 1;
            rc = -5;
          end
        end
      else % hpdeg < degree
        [zFar,zFarErr,zUC,zUCErr,sp_zFar,sp_zUC,LNew]=...
                         updateR(c,cErr,zFar,zFarErr,zUC,zUCErr,sp_zFar,sp_zUC);
        if(~LNew) % c contained no new roots.
          z = [zFar ; zUC];
          zErr = [zFarErr; zUCErr];
          done = 1;
          rc = -6;
        end
      end
    end
  end
end

% I did tlroots(1000,1:60000) and 16 failed with rc=-2, "correct degree, but
% unfactor() did not confirm." It did a second grid search and then it worked. 
% That was unnecessary. It found all roots the first time! I looked at the first
% in detail. Array c had 6 complex roots with arguments 105 145 168 172 173.07 173.09 degrees.
% Four were close together and close to the x axis so they were close to their
% conjugates. Three were extremely close so the dynamic range of sp_zFar was 3e+11.
% It was enormous near 173 degrees because it was missing 3 roots that would have
% suppressed the spectrum there. maxErr above was about 5e-4. I tried computing u above
% two other ways: 1) compute u completely in the time domain 2) compute spectrum with
% fft(poly(leja(hp2fp(zFar))),1024) etc. The resulting u's still had polycmp error
% of about 5e-4. However, I also did 
% z=[zFar;zUc;c]; u=unfactor(z,1); v=poly(leja(hp2fp(z))); Then polycmp(x,u)=1e-13 and 
% polycmp(x,v)=1e-14. Using leja was 9 times slower. It seems I could only get an accurate
% unfactorization if I first combined all the roots. BTW, zUC was 1.0. After adding
% the following section, all 16 worked the first time.

if(rc == -2)
   u = unfactor(z,1);
   maxErr = polycmp(x,u);
   if(maxErr < max_error_allowed)
      rc = nDeflations;
%      hwrite('** FYI I needed to use the "rc == -2" fixup in iterativeDeflation().');
   end
end

% BATCH MODE: disabled diagnostic output
% if(~isempty(deflatedR))
%   fp = fopen(DEFLFILE, 'a');
%   if(fp)
%     fprintf(fp, '\n%g %g %g %g %g %g %g %g %g %g', sort(deflatedR));
%     fclose(fp);
%   end
% end

if(nargout == 4)
  outParms = [length(deflatedR), maxErr];
end

if(NZeros)
   z    = [z; zeros(NZeros,1)];
   zErr = [zErr; zeros(NZeros,1)];
end

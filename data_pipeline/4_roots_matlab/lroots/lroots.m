% function [rts, rtsErr, rc, outParms] = lroots(x, inParms)
%
% PURPOSE:
%   Polynomial factorization by the Lindsey-Fox grid-seach-in-polar-coordinates method.
%
% REQUIRED INPUT
%   x - the polynomial coefficients.
%
% OPTIONAL INPUT
%   inParms - Optional input parameters. inParms can have anywhere from 1
%             to 7 values. If you omit a value, you must omit all the
%             following values.
%     inParms(1) - ID number for information messages in the history file. 
%                  0 suppresses messages. 'Messages' means both error and 
%                  success messages. Default=1.
%     inParms(2) - Default=1. A value > 1 determining radial extent. Larger
%                  values extend the region with a fine grid size to radii
%                  farther from 1. For example, Butterworth filters have roots
%                  far from the unit circle. One could be factored with lroots(b,[1,2]);
%     inParms(3) - Default=1. A power of 2 determining angular grid size.
%                  Larger values decrease the grid size.
%     inParms(4) - Default=1. A value > 1 determining radial grid size.
%                  Larger values decrease the radial grid size.
%     inParms(5) - Default=5.  The first break where we double the angular grid size.
%     inParms(6) - Default=10. The second break where we double the angular grid size.
%     inParms(7) - Default=18. The third break where we double the angular grid size.
%
% OUTPUT
%   rts - The output array of roots as a column vector.
%   rtsErr - Error estimates for the roots.
%   rc - The return code. See iterativeDeflation for the return codes.
%   outParms - An array of output parameters if the deflation succeeded.
%              outParms(1) = the number of roots obtained by deflation.
%              outParms(2) = max(abs(|x - unfactor(z)|)).

% USES
%   extremeRoots, gridSearch, hpdegree, iterativeDeflation, ixFirstNZ, ixLastNZ,
%   laguerre3, lastValue, nonPolish, normalizePow2, nroots, quadraticRoots, radiusGen, uniq

% AUTHOR
%   James W. Fox 4/24/2001

% COMMENT
% -------
% The FFT can be used to efficiently evaluate any n'th degree polynomial at 
% r*(the n'th roots of unity) where r is any real scalar. Joe Pat Lindsey
% observed that this allows an efficient grid seach in polar coordinates.

% James W. Fox implemented his ideas and added several ideas about polishing
% roots with very large magnitude, accurately unfactoring roots, and recursively
% deflating roots.

function [rts, rtsErr, rc, outParms] = lroots(x, inParms)

tic;

% Declare globals
global DIAGNOSTICS RAM MAX_ARRAY_LENGTH_EXPONENT HISTORYFILE DEFLFILE;

% Initialize globals
rootsIni;

s = size(x);
if(s(1) ~=1 & s(2) ~= 1)
   error('$1 must be 1xn or nx1.')      
end

ver = str2double(regexp(version,'^\d+','match','once'));

%---- Error checks ----
if(isempty(x))
  rts = [];
  rc = 0;
  outParms = [0,0];
  return;
elseif(ischar(x))
   error('$1 must be a numeric array.');
elseif(max(isnan(x)))
   error('$1 must not contain NaN or -NaN.');
elseif(max(isinf(x)))
   error('$1 must not contain Inf or -Inf.');
elseif(length(x)==1 & x~=0)
  rts = [];
  rc = 0;
  outParms = [0,0];
  return;
elseif(~max(abs(x)))
  % BATCH MODE: disabled diagnostic output
  % fprintf('Polynomial is entirely zero. Every number is a root.\n');
  return;
end

if(nargin == 1) 
  inParms = 1;
end
[ID_NUM, radialExtent, thetaMult, radialMult, WindowBreak1, WindowBreak2, WindowBreak3] = processInparms(inParms);

%---- Remove leading and trailing zeros and normalize ----
ixLNZ = ixLastNZ(x);
nZeros = length(x) - ixLNZ; % Multiplicity of 0.0 as a root. Used at the very end.
x = normalizePow2(x(ixFirstNZ(x) : ixLNZ));
Lx = length(x);
degree = Lx - 1; % Note, degree == 0 is possible. Consider x=[0 0 1 0]

if(degree <= 2)
   if(degree == 0) % If originally x=[0 0 1 0] then the above changes it to x=1
      rts = [];
   elseif(degree == 1)
      rts = -x(2)/x(1);
      rtsErr = 0;
   elseif(degree == 2)
      [rts,rtsErr] = laguerre3(x, quadraticRoots(x));
   end
   
   if(nZeros)
      rts = [rts; zeros(nZeros,1)];
      rtsErr = [rtsErr; zeros(nZeros,1)];   
   end
   rc = 0;
   outParms = [0,0];
   return;
end

exponGain = 1;      % Necessary if I comment out expGain().
xx = x;             % Necessary if I comment out expGain().
extremeRts = [];    % Necessary if I comment out extremeRoots()
extremeRtsErr = []; % Necessary if I comment out extremeRoots()

% I have seen exponential gain sabotage the search for extreme roots which are sensitive.
% Therefore do not do exponential gain before extremeRoots as I originally did.
%
% Well, there are two sides to this. If the polynomial is exponentially increasing and
% has huge coefficients and a huge root, then polynomial evaluation may overflow. However,
% if you first do exponential gain, then it would be less likely to overflow. I am not
% sure which is more important.
[extremeRts, extremeRtsErr, xx, rc, outParms] = extremeRoots(x);
if(rc >= 0)
   rts = extremeRts;
   rtsErr = extremeRtsErr;
   % We previously stripped off trailing zeros.
   if(nZeros)
      z = [rts; zeros(nZeros,1)];
      zErr = [rtsErr; zeros(nZeros,1)];
   end
   [rts, rtsErr] = hp2fp(rts, rtsErr);
   return
end

[xx, exponGain] = expGain(xx(:).');

% The above line finds AND deflates extreme roots. If some of them are ill-conditioned,
% the resulting polynomial may not be accurate. Perhaps it would be best to find these
% roots without deflating them. Just add them to roots found by grid search. However,
% be careful. Suppose 1e-25*(4th roots of unity) are extreme roots. Be careful that
% you do not declare them to be false double roots. The current code will do that because
% you do not use the accuracy of each root when determining double roots.
xxflipped = flip(xx);

% Note, extremeRts are roots of x but not xx. rts are roots of xx.

%---- Determine the final radial and theta multipliers ----
degree = length(xx) - 1;
[radialM, thetaM] = findMultipliers(degree, radialMult, thetaMult);
inc = 1/degree;

%------ Start Searches at 0 and at infinity ------

% Grid search will not find roots with tiny radius.
% Fill in the gap by polishing 0.0 and seeing where it lands.
% Since we start at (0,0) we might polish to a root in the lower 
% half plane. Ensure we only get upper half plane roots.

% Note, if the linear and quadratic terms are 0 then f'(0)=f''(0)=0
% and the denominator of the laguerre correction is 0. I have not done
% anything in the laguerre code to deal with this. Therefore, I am going
% to start at 1e-100*i instead of starting at 0.

rts = [];

rMax=[];
[z, zErr] = laguerre3(xx, 1e-100*i);
z = nonPolish(z, zErr);
if(length(z))
  rMax = z;
  if(imag(rMax)<0) % We only want upper half-plane roots.
    rMax = conj(rMax);
  end
end

rMin=[];
[z, zErr] = laguerre3(xxflipped, 1e-100*i);
z = nonPolish(z, zErr);
if(length(z))
  rMin = 1/z;   % nonPolish() ensures z~=0 or Inf.
  if(imag(rMin)<0) % We only want upper half-plane roots.
    rMin = conj(rMin);
  end
end

% Someday in the future, if you find very small or very large roots,
% repeat the above search to see if there are any more. Then rMin,
% and rMax may have more than 1 value. Eg.
% while(10*abs(x(1)) < abs(x(2))) % continue looking for large magnitude roots.

rts = [rMin;rMax];

% By the way, rMin and rMax will probably not be empty. If there is a lone
% extreme root it will be found. If not, they will find some root with
% normal magnitude which will probably be found a second time by grid
% search. I will rely on uniq() to get rid of the second copy.

%------------------ Use Grid Search -------------------

% You need a better criterion for when to give up. Currently if I factor
% a 1000 degree poly and repeat a root to make it a root of order 4,
% the following loops until it reaches the maximum.

rc = -Inf; % Signals this is the first time through the while loop.
NPasses = 0;   % The number of grid searches already performed.
saveTime = toc;

while(rc < 0 & rc ~= -1 & NPasses < 4 & thetaM*2^nextpow2(length(xx)) <= 2^MAX_ARRAY_LENGTH_EXPONENT)
  NPasses = NPasses + 1;
  if(rc ~= -Inf) % If this is not the first time, decrease grid size.
    % For random coefficient polys, if I fail on the first try, doing 
    % 'thetaM = 2*thetaM' once has always worked on
    % the second pass. I did not need to also increase radialM. Thus I 
    % only increase radialM if we get to 3 passes. This will happen if 
    % you try to factor a Butterworth using defaults. Butterworths
    % require 3 or 4 passes if you use defaults.
   
    thetaM = 2*thetaM;
    if(NPasses >= 3)
      radialExtent = 2*radialExtent;
      radialM = 2*radialM;
    end
   
    if(DIAGNOSTICS)
      % BATCH MODE: disabled diagnostic output
      % fprintf('\nthetaM=%d \tradialM=%d \tRAM=%d\n',thetaM, radialM, RAM);
      % hwrite(' ');
      % hwrite(['thetaM=',num2str(thetaM),', radialM=',num2str(radialM),', RAM=',num2str(RAM)]);
    end
    
    if(degree >= 2^18 & ver >= 6) % pack can be called in version 6 or 7 but not in 5.
      packMemory;
    end
  end

  %----------------------------------------
  searchRadii = radiusGen(xx, radialExtent, radialM); % An efficient set of search radii.
  ixWindowBreak1 = lastValue(find(searchRadii >= 1-WindowBreak1*inc*radialExtent));
  ixWindowBreak2 = lastValue(find(searchRadii >= 1-WindowBreak2*inc*radialExtent));
  ixWindowBreak3 = lastValue(find(searchRadii >= 1-WindowBreak3*inc*radialExtent));
  
  rts1 = zeros(degree, 1);
  nRts1 = 0; % The number of values saved in rts1
  %----------------------------------------
  % This section checks for roots at r==1 by checking at only 3 radii.
  % xxflipped is not used.
  radii = [searchRadii(2:-1:1), 1/searchRadii(2)];
  thisTheta = max(1, thetaM); 
  if(length(radii) >= 3)
    [rts1, nRts1] = gridSearch(xx, radii, thisTheta, 0, rts1, nRts1); % Roots near Unit Circle
  end
  %----------------------------------------
  radii = searchRadii(ixWindowBreak1 : -1 :  1);
  if(length(radii) >= 3)
    [rts1, nRts1] = gridSearch(xxflipped, radii, thisTheta, 1, rts1, nRts1);
    [rts1, nRts1] = gridSearch(xx, radii, thisTheta, 0, rts1, nRts1);
  end
   %----------------------------------------
  radii = searchRadii(ixWindowBreak2 : -1 :  ixWindowBreak1-1);
  thisTheta = max(1, thetaM/2);
  if(length(radii) >= 3)
    [rts1, nRts1] = gridSearch(xxflipped, radii, thisTheta, 1, rts1, nRts1);
    [rts1, nRts1] = gridSearch(xx, radii, thisTheta, 0, rts1, nRts1);
  end
  %----------------------------------------
  radii = searchRadii(ixWindowBreak3 : -1 : ixWindowBreak2-1);
  thisTheta = max(1, thetaM/4);

  if(length(radii) >= 3)
    [rts1, nRts1] = gridSearch(xxflipped, radii, thisTheta, 1, rts1, nRts1);
    [rts1, nRts1] = gridSearch(xx, radii, thisTheta, 0, rts1, nRts1);
  end
  %----------------------------------------
  radii = searchRadii(length(searchRadii) : -1 : ixWindowBreak3-1);
  thisTheta = max(1, thetaM/8);

  if(length(radii) >= 3)
    [rts1, nRts1] = gridSearch(xxflipped, radii, thisTheta, 1, rts1, nRts1);
    [rts1, nRts1] = gridSearch(xx, radii, thisTheta, 0, rts1, nRts1);
  end
  rts1 = rts1(1:nRts1); % Truncate to how many we actually found.
  
  if(DIAGNOSTICS)
    % BATCH MODE: disabled diagnostic output
    % fprintf('<Grid searching polynomial # %d found %d putative upper half-plane roots and took %g seconds>\n', ID_NUM, length(rts1), toc);
    % hwrite(['<Grid searching polynomial #', num2str(ID_NUM), ' found ', num2str(length(rts1)), ...
    %          ' upper half-plane roots and took ', num2str(toc-saveTime), ' seconds.']);
  end

  %----------------------------------------
  
% rts may equal [], but it may contain roots if 1) we obtained rts from 
% iterativeDeflation above or 2) this is the second pass through the loop.
  rts = [rts; rts1]; 
  
  saveTime = toc;
  [rts, rtsErr] = bigPolish(xx, rts);
% [rts, rtsErr] = bigPolish(x, rts);  % 2/22/2006 This can find extreme
% roots which are not roots of xx so iterativeDeflation screws up.
  
  [rts, rtsErr] = nonPolish(rts, rtsErr);
  expF = 1; % Fudge factor. Use 12 if you use relCorr error estimates; they are smaller.
  [rts, rtsErr] = uniq(rts, expF*rtsErr);
  rtsErr = rtsErr/expF;
  % If rtsErr is ever large and expF>>1, then uniq may eliminate several legitimate roots.
  
  NLeft = degree - hpdegree(rts); % The exact (if uniq worked) number of missing roots. 
  if(NLeft < 0) % Oops. uniq() failed. Try polishing again. This worked once.
      
%  [rts, rtsErr] = laguerre3(x, rts); % 2/22/2006 See comment above.
   [rts, rtsErr] = laguerre3(xx, rts);

    [rts, rtsErr] = nonPolish(rts, rtsErr); % the starting position was totally off.
    [rts, rtsErr] = uniq(rts, []);
    NLeft = degree-hpdegree(rts);
  end
  if(DIAGNOSTICS)
    % BATCH MODE: disabled diagnostic output
    % fprintf('<Polishing and removing duplicates left %d roots missing and took %g seconds>\n', NLeft, toc-saveTime);
    % hwrite(['<Polishing and removing duplicates left ', num2str(NLeft), ...
    %         ' roots missing and took ', num2str(toc-saveTime), ' seconds.']);
  end

  % The folowing added by Ryan King to help keep track of progress on large polynomials.
  if (degree >= 2^17)
     % BATCH MODE: disabled diagnostic output
     % fprintf('%s\n',datestr(now))
  end;
  
  ver = str2double(regexp(version,'^\d+','match','once'));
  if(degree >= 2^18 & ver >= 6) % pack can be called in version 6 or 7 but not in 5.
    % BATCH MODE: disabled diagnostic output
    % fprintf('Defragmenting memory before deflation...   ');
    packMemory;
    % fprintf(' done.\n');
  end
    
  expF = 1; % Ad hoc expansion factor. This should be small.
  saveTime = toc;
  [rts, rtsErr, rc, outParms] = iterativeDeflation(xx, rts, expF*rtsErr, ID_NUM);
  rtsErr = rtsErr/expF;
  if(DIAGNOSTICS)
    % BATCH MODE: disabled diagnostic output
    % fprintf('<Unfactoring and deflating to find missing roots took %g seconds>\n', toc-saveTime);
    % hwrite(['<Unfactoring and deflating to find missing roots took ', num2str(toc-saveTime), ' seconds.']);
  end

  if(ID_NUM)   % If ID_NUM is non-zero, print a success or failure message.
    % BATCH MODE: disabled diagnostic output
    % if(rc == 0)   % No deflation was necessary.
    %    hwrite(['Already found all roots for poly #', num2str(ID_NUM), ', seconds from start=',num2str(toc), ',  Max error=',num2str(outParms(2))]);
    % elseif(rc == -1)
    %   hwrite(['** Failure for poly #', num2str(ID_NUM), ...
    %           ', too many roots. Probable polish error.']);
    %   hwrite('** ==============> Impossible to continue <==============');
    % elseif(rc == -2)
    %   hwrite(['** Failure for poly #', num2str(ID_NUM), ...
    %           ' correct degree, but unfactor() did not confirm.']);
    % elseif(rc == -4)
    %    hwrite(['** Failure for poly #', num2str(ID_NUM), ...
    %          ' Deflation provided no new roots.']);
    % elseif(rc == -5)
    %   hwrite(['** Failure for poly #', num2str(ID_NUM), ...
    %          ' correct degree, but unfactor() did not confirm and no new roots.']);
    % elseif(rc == -6)
    %   hwrite(['** Failure for poly #', num2str(ID_NUM), ...
    %          ' degree too small and no new roots.']);
    % else       % Deflation worked
    %   hwrite(['<Poly #', num2str(ID_NUM), ': ', num2str(rc), ...
    %         ' deflation(s) for ',num2str(outParms(1)), ' upper half-plane root(s). Seconds from start=', num2str(toc), ...
    %         '. Max error=', num2str(outParms(2)),'>']);
    %   if(DIAGNOSTICS)
    %     hwrite(' ');
    %   end   
    % end
  end

  % Normally the next line should be commented out. Uncomment it to keep it from
  % doing retries upon deflation failure. I.e., this allows you to analyze 
  % deflation failures.
  % This is dangerous!!! If you uncomment it, you may fail to
  % put the comment back. You would think every factorization
  % was working even though sometimes it was not.
  %          rc = 0; % Normally this should be commented out.
end

% fprintf('Poly # %d: %d deflation(s) for %d upper half-plane roots. Seconds =%g. Error=%g.\n', ID_NUM, rc, outParms(1), toc, outParms(2));

%-------------------------------------------------------
[junk,j] = sort(angle(rts));
rts = exponGain*rts(j);
rtsErr = rtsErr(j);

% extremeRts were roots of x but not xx so we can not append them while working with xx above.
rts = [rts; extremeRts];
rtsErr = [rtsErr; extremeRtsErr];

% We previously stripped off trailing zeros.
if(nZeros)
   rts = [rts; zeros(nZeros,1)];
   rtsErr = [rtsErr; zeros(nZeros,1)];
end

% Up to this point we have only computed upper-half plane roots.
% Now supply lower half-plane roots too.
%
% FIXED: snap near-real roots to exactly real first, using the same
% tolerance convention as the hpdegree.m fix above. A numerically-
% polished real root almost never comes back with imag() exactly 0.0,
% so the strict imag(rts)>0 test below would otherwise treat it as
% complex and mirror a spurious near-duplicate (its own tiny-residue
% conjugate) into the lower half plane instead of leaving it as the
% single real root it actually is.
REALTOL = 1e-8;
isNearReal = abs(imag(rts)) <= REALTOL * max(1, abs(rts));
rts(isNearReal) = real(rts(isNearReal));

k=find(imag(rts) > 0);
rts = [rts; conj(rts(k))];
rtsErr = [rtsErr; rtsErr(k)];

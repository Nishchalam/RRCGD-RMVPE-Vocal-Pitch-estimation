% function [z, x, zErr, RC, parms] = tlroots(degree, seeds, inParms);
%
% PURPOSE:
%   Driver to repeatedly test lroots().
%
% REQUIRED INPUT:
%   degree - The degree of the random coefficient polynomials to be generated.
%   seeds - An array of integers to be passed as seeds to rand() and randn().
%
% OPTIONAL INPUT:
%   inParms - Optional third argument to lroots.
%
% OUTPUT:
%   z - The roots of x.
%   x - The polynomial.
%   Information to files roots.his and defl.his

% USES:
%   hwrite, rootsIni, lroots
%
% AUTHOR:
%   James W. Fox 4/25/2001


function [z, x, zErr, RC, parms] = tlroots(degree, seeds, inParms);

global DIAGNOSTICS DEFLFILE;

% Initialize globals.
rootsIni;

% Erase the list of deflated radii.
fid=fopen(DEFLFILE);
if(fid ~= -1) % If it exists.
   fclose(fid);
   delete(DEFLFILE);
end

if(nargin < 2)
  error('Supply $2.');
elseif(length(degree)~=1 | degree ~= round(degree) | degree < 1)
  error('$1 must be a postive integer.')
elseif(any(seeds - round(seeds) | seeds < 0))
   error('$2 must only contain non-negative integers.')
elseif(nargin < 3)
   inParms = 1;
end

hwrite(' ');
if(length(seeds) <= 10)
   if(nargin == 2)
     hwrite(['   tlroots(',num2str(degree),',[',num2str(seeds),'])'])
   else   
     hwrite(['   tlroots(',num2str(degree),',[',num2str(seeds),'],[',num2str(inParms),']'])
   end
else
   if(nargin == 2)
     hwrite(['   tlroots degree=',num2str(degree),' and # of seeds = ',num2str(length(seeds))]);
   else
     hwrite(['   tlroots degree=',num2str(degree),' and # of seeds = ',num2str(length(seeds)),' and inParms=[',num2str(inParms),']']);
   end
end

% rc(k) and NDeflated(k) will ultimately contain the return code and
% the number of deflated roots for the polynomial using seeds(k).

Lseeds = length(seeds);
rc = zeros(1,Lseeds);
NDeflated = zeros(1,Lseeds);
unfactorErr = zeros(1,Lseeds);
startTime = clock;
for k=1:Lseeds
%  fprintf('\n\t<<<<< tlroots beginning seed = %d >>>>>\n', seeds(k)) ;
  rand('seed', seeds(k));
  x = rand(1, degree+1);
  x = normalizeJWF(x - mean(x)); % uniformly distributed between -.5 and .5.
  %  x = x + .001; % To make sure mean was no 0 exactly. This is no longer important.
  inParms(1) = seeds(k);
  [z, zErr, RC, parms] = lroots(x, inParms);
  rc(k) = RC;
  unfactorErr(k) = parms(2);
  NDeflated(k) = parms(1);
end

seconds = etime(clock,startTime);
k = find(rc >= 0);
if(length(k)) 
 NDeflated = NDeflated(k);
end

NFailures = length(rc) - length(k);

if(~DIAGNOSTICS)
  hwrite(' ');
end
hwrite(['Average seconds including overhead = ',num2str(seconds/Lseeds)]);
hwrite(['Max number of deflated upper half plane roots = ',num2str(max(NDeflated))]);
hwrite(['Average number of deflated upper half plane roots = ',num2str(mean(NDeflated))]);
hwrite(['Average unfactor error = ',num2str(mean(unfactorErr))]);
if(NFailures)
  hwrite(['Number of total failures = ',num2str(NFailures),' <===============']);
end


fprintf('\nAverage seconds including overhead = %g\n',seconds/Lseeds);
fprintf('Max number of deflated upper half plane  roots = %g\n',max(NDeflated));
fprintf('Average number of deflated upper half plane roots = %g\n',mean(NDeflated));
fprintf('Average unfactor error = %g\n',mean(unfactorErr));
if(NFailures)
fprintf('Number of total failures = %g<===============\n',NFailures);
end

% function [polished, x, raw, err, rc, parms] = tbRoots(degree, seeds);
%
% PURPOSE
%   Driver to repeatedly test bRoots().
%
% REQUIRED INPUT
%   degree - The degree of the random coefficient polynomials to be generated.
%   seeds - An array of integers to be passed as seeds to rand().
%
% OUTPUT
%   x        - The last polynomial
%   raw      - Unpolished roots of the last polynomial
%    polished - Polished roots of the last polynomial
%   err      - Error estimates for the roots.
%   RC       - An array of return codes.
%   parms    - DeflateFFT parms for the last polynomial.
%   Information written to files roots.his and defl.his

% USES
%   rootsIni, normalizeJWF
%
% AUTHOR
%   James W. Fox 10/9/2002

function [polished, x, raw, err, rc, parms] = tbRoots(degree, seeds);

global HISTORYFILE DEFLFILE;

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
elseif(degree ~= round(degree) | degree < 1)
  error('$1 must be a postive integer.')
elseif(any(seeds - round(seeds) | seeds < 0))
  error('$2 must only contain non-negative integers.')
end

hwrite(' ');
hwrite([' 		broots: degree=', num2str(degree)]);

% rc(k) and NDeflated(k) will ultimately contain the return code and
% the number of deflated roots for the polynomial using seeds(k).

Lseeds = length(seeds);
rc = zeros(1,Lseeds);
NDeflated = rc;
startTime = clock;
k = 1;
while(k <= Lseeds)
%  fprintf('\n\t<<<<< tbroots beginning Seed = %d >>>>>\n', seeds(k)) ;
  rand('seed',seeds(k));
  x = rand(1,degree+1);
  x = normalizeJWF(x - mean(x)); % uniformly distributed between -.5 and .5.
  x = x + .001; % To test the setting of thetaK = 1 or 2.
  [raw, polished, err, RC, parms] = broots(x, seeds(k));
  rc(k) = RC;
  NDeflated(k) = parms(1);
  k = k+1;
end

seconds = etime(clock,startTime);
k = find(rc >= 0);
if(length(k)) 
 NDeflated = NDeflated(k);
end

NFailures = length(rc) - length(k);
if(NFailures)
  hwrite(['Number of total failures = ',num2str(NFailures),' <===============']);
end
hwrite(['Average seconds = ',num2str(seconds/Lseeds)]);
hwrite(['Max number of deflated roots = ',num2str(max(NDeflated))]);
hwrite(['Average number of deflated roots = ',num2str(mean(NDeflated))]);

if(NFailures)
fprintf('\n%s',['Number of total failures = ',num2str(NFailures),' <===============']);
end
fprintf('\n%s',['Average seconds = ',num2str(seconds/Lseeds)]);
fprintf('\n%s',['Max number of deflated roots = ',num2str(max(NDeflated))]);
fprintf('\n%s\n',['Average number of deflated roots = ',num2str(mean(NDeflated))]);

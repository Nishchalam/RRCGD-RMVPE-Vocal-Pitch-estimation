% function [err1, err2,x, z, zErr, RC, parms] = tunrap(degree, seeds);
%
% PURPOSE
%   Compare unwrapm and unwrapph with unwrapz (which is the correct unwrapping).
%
% REQUIRED INPUT
%   degree - The degree of the random coefficient polynomials to be generated.
%   seeds - An array of integers to be passed as seeds to rand() and randn().
%
% REQUIRED OUTPUT
%   err1 - The maximum absolute error between unwrapz and unwrapph.
%
% OPTIONAL OUTPUT
%   err2 - The maximum absolute error between unwrapz and unwrapm.

% USES:
%   rootsIni, lroots, unwrapm, unwrapph, unwrapz
%
% AUTHOR:
%   James W. Fox 1/3/2004


function [err1, err2, x, z, zErr, RC, parms] = tunwrap(degree, seeds);

global HISTORYFILE DEFLFILE;

% Initialize globals.
rootsIni;

% Erase the list of deflated radii.
delete(DEFLFILE);

if(nargin < 2)
  error('Supply $2.');
elseif(degree ~= round(degree) | degree < 1)
  error('$1 must be a postive integer.')
elseif(any(seeds - round(seeds) | seeds < 0))
  error('$2 must only contain non-negative integers.')
end

hwrite(' ');
hwrite([' 	lroots: degree=', num2str(degree)]);

% rc(k) and NDeflated(k) will ultimately contain the return code and
% the number of deflated roots for the polynomial using seeds(k).

Lseeds = length(seeds);
err1 = zeros(1,Lseeds);
err2 = err1;
rc = err1;
NDeflated = err1;
unfactorErr = err2;
startTime = clock;
for k=1:Lseeds
%  fprintf('\n\t<<<<< tunwrap beginning seed = %d >>>>>\n', seeds(k)) ;
  rand('seed', seeds(k));
  x = rand(1, degree+1);
  x = x - mean(x); % Ensure mean==0.
  x = x + .001; % If mean==0, only unwrapz will get the correct answer.
  [z, zErr, RC, parms] = lroots(x, seeds(k));
  rc(k) = RC;
  unfactorErr = parms(2);
  NDeflated(k) = parms(1);
  p = unwrapz(x,z,zErr);
  ph = unwrapph(x);
  err1(k) = max(abs(p-ph));
  if(nargout > 1)
     m = unwrapm(x);
     err2(k) = max(abs(p-m));
  end
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
hwrite(['Max number of deflated upper half plane roots = ',num2str(max(NDeflated))]);
hwrite(['Average number of deflated upper half plane  roots = ',num2str(mean(NDeflated))]);
hwrite(['Average unfactor error = ',num2str(mean(unfactorErr))]);

if(NFailures)
fprintf('Number of total failures = %g<===============\n',NFailures);
end
fprintf('Average seconds = %g\n',seconds/Lseeds);
fprintf('Max number of deflated upper half plane  roots = %g\n',max(NDeflated));
fprintf('Average number of deflated upper half plane  roots = %g\n',mean(NDeflated));
fprintf('Average unfactor error = %g\n',mean(unfactorErr));

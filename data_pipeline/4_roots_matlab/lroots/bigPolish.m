% functiono [rts, rtsErr] = bigPolish(xx, rts);
%
% PURPOSE:
%   Polish a large number of roots.
%
% REQUIRED INPUT:
%   xx - the polynomial coefficients
%   rts - the starting positions.
%
% OUTPUT:
%   rts - the polished positions
%   rtsErr - error estimates for the roots.

% USES:
%   laguerre3
%
% AUTHOR:
%   James W. Fox 1/31/2006

function [rts, rtsErr] = bigPolish(xx, rts);

global nCPUs;

if(nCPUs == 1) % There is only 1 cpu.
  [rts, rtsErr] = laguerre3(xx, rts);
else
   % Partition rts into the first 1/nCPUs roots, and the second 1/nCPUs roots etc. Then send them to the nCPUs 
   % different CPUs and call laguerre3 nCPUs times in parallel. Concatenate the results into rts and rtsErr.
   %
   % bucketSize=max(1, floor(length(rts)/nCPUs));
   % for k=0:nCPUs-2 % All but the last bucket
   %   rtsPartial = rts(k*bucketSize+1:(k+1)*bucketSize);
   %   Send xx and rtsPartial to another cpu and call laguerre3.
   % end
   %
   % rtspartial = rts((nCPUs-1)*bucketSize+1:end); % The last bucket may be an odd length:
   % Send xx and rtsPartial to another cpu and call laguerre3.
end

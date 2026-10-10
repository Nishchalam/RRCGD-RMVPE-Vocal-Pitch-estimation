% function [zFar,zFarErr,zUC,zUCErr,sp_zFar,sp_zUC,LNew] = ...
%          updateR(c,cErr,zFar,zFarErr,zUC,zUCErr,sp_zFar,sp_zUC)
%
% PURPOSE
%   Suppose we already have an incomplete set of polished roots, zFar & zUC.
%   After doing polynomial division, we MAY have found new roots, c.
%   However, c may contain roots previously found in zFar. If c contains
%   any new roots, append them to zFar and zUC and update their spectra. 
%   Also return the number of new roots.
%
% REQUIRED INPUT:
%   c - A tentative set of new roots. Duplicates will be eliminated.
%   cErr - An array of error estimates for c. If only a single number is supplied,
%          it is repeated. If [] is supplied, 1e-7 is used.
%   zFar - A previously generated incomplete set of polished roots not near the unit circle.
%   zFarErr - An array or error estimates for zFar. If only a single number is supplied,
%          it is repeated. If [] is supplied, 1e-7 is used.
%   zUC - A previously generated incomplete set of polished roots near the unit circle.
%   zUCErr - An array of error estimates for zUC. If only a single number is supplied,
%          it is repeated. If [] is supplied, 1e-7 is used.
%   sp_zFar - the spectrum of unfactor(zFar,1).
%   sp_zUC - the spectrum of unfactor(zUC,1).

%   Note, neither zFar nor zUC will be checked for duplicates.
%   Note, if either zFar or zUC is empty, the spectrum must be all 1's and not [].
%
% OUTPUT:
%   Updated values for the last 6 input arguments and
%   LNew - the number of new values added from c. It may be 0.

% USES:
%   ismembr, uniq
%
% AUTHOR:
%   James W. Fox 5/8/2003

function [zFar,zFarErr,zUC,zUCErr,sp_zFar,sp_zUC,LNew] = ...
   updateR(c,cErr,zFar,zFarErr,zUC,zUCErr,sp_zFar,sp_zUC)

if(isempty(cErr))
   cErr = 1e-7 + zeros(size(c));
elseif(length(cErr) == 1)
   cErr = abs(cErr) + zeros(size(c));
else                          
  if(length(cErr) ~= length(c))
    error('length($1) must equal length($2).')
  elseif (min(cErr) < 0)
    error('$2 must be not contain negative values.')
  end
end

if(isempty(zFarErr))
   zFarErr = 1e-7 + zeros(size(zFar));
elseif(length(zFarErr) == 1)
   zFarErr = abs(zFarErr) + zeros(size(zFar));
else                          
  if(length(zFarErr) ~= length(zFar))
    error('length($3) must equal length($4).')
  elseif (min(zFarErr) < 0)
    error('$4 must be not contain negative values.')
  end
end

if(isempty(zUCErr))
   zUCErr = 1e-7 + zeros(size(zUC));
elseif(length(zUCErr) == 1)
   zUCErr = abs(zUCErr) + zeros(size(zUC));
else                          
  if(length(zUCErr) ~= length(zUC))
    error('length($5) must equal length($6).')
  elseif (min(zUCErr) < 0)
    error('$6 must be not contain negative values.')
  end
end

Lsp = length(sp_zFar); 
[c,cErr] = uniq(c(:), cErr); % Root arrays, like zFar, are column vectors. Ensure c is also.
j = find(~ismembr(c,cErr,[zFar;zUC],[zFarErr;zUCErr])); % Indicies of new roots. May be empty.
c = c(j); % New roots not found in zFar or zUC. Might be empty.
cErr = cErr(j);
LNew = length(j);
if(~LNew) % No new roots to add.
  return
end

test  = abs(abs(c)-1) < 1e-12; % True/false where distance to unit circle is small/large.
ixUC  = find(test);  % Indicies of roots close to the unit circle.
ixFar = find(~test); % Indicies of roots far from the unit circle.
cUC   = c(ixUC);     % New roots close to the unit circle.
cFar  = c(ixFar);    % New roots far from the unit circle.

if(~isempty(ixUC))
  zUC = [zUC ; cUC];
  zUCErr = [zUCErr; cErr(ixUC)];
  [junk, sp_UC]= unfactor(cUC, 1, 0, Lsp);
  sp_zUC = sp_zUC .* sp_UC;
end

if(~isempty(ixFar))
  zFar = [zFar ; cFar];
  zFarErr = [zFarErr; cErr(ixFar)];
  [junk, sp_cFar] = unfactor(cFar, 1, 0, Lsp);
  sp_zFar = sp_zFar .* sp_cFar;
end

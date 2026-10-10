% function [rts1, nRts1] = gridSearch(c, radii, padfactor, invert, rts1, nRts1);
%
% PURPOSE
%   Use the FFT to do a fast grid search to find approximate root locations.
%
% REQUIRED INPUT:
%   c - The polynomial coefficients. c(1)=the constant term.
%       c must not have leading or trailing zeros. We may want to
%       invert the radii and if so, 0 must not be a root.
%   radii - An array of search radii.  length(radii) >=  3.
%   padfactor - A small multiple of 2. Increasing this value
%               makes the search grid finer in the theta direction.
%   invert - If 1, invert the roots. If 0, do not invert. The driver
%            searches for roots outside the unit circle by searching
%            for roots of the flipped polynomial that are iside the
%            unit circle. After they are found, and polished, they
%            must be inverted.
%
% OUTPUT:
%   rts1 - A column vector of new polished roots found by searching
%           this radial range.
%   nRts1 - The number of values stored in rts1.

% USES:
%   newton3, normalizePow2, nonPolish, rdft
%
% COMMENT:
%   This function creates 4 long arrays, cPadded, F1, F2, and F3.
%
% AUTHOR:
%   James W. Fox 4/22/2001
%
% CHANGES
%   12/19/2005 Write FFT's out to 3 files for extremely large polynomials. 
%   1/9/2006   Add call to packMemory

function [rts1, nRts1] = gridSearch(c, radii, padfactor, invert, rts1, nRts1);

global DIAGNOSTICS
if(isempty(DIAGNOSTICS))
  DIAGNOSTICS = 0;
end

DEGREE_FOR_USING_FILES = 2^19; % Save FFT's to 3 files if degree >= this.
degree = length(c) - 1;
Nradii = length(radii);

ver = str2double(regexp(version,'^\d+','match','once'));
if(degree >= 2^20 & ver >= 6) % pack can be called in version 6 or 7 but not in 5.
  clear ver;
  packMemory;
end

c = c(:).'; % Necessary to avoid error msg. "Matrix dimensions must agree."
if(~c(1));
  error('$1(1) must not be zero.');
elseif(~c(length(c)));
  error('The last value of $1 must not be zero.');
elseif(degree == 0);
  error('$1 must have at least 2 coefficients.');
elseif(Nradii < 3);
  error('$2 must have at least 3 search radii.');
elseif(padfactor < 1);
  error('$3 must be >= 1.');
end

if(degree >= DEGREE_FOR_USING_FILES)
   fpName1 = tempname;
   fpName2 = tempname;
   fpName3 = tempname;
   
   fp = zeros(1,3);
   fp(1) = fopen(fpName1,'w+');
   if(fp(1) == -1)
      error('Could not open first temp file to hold FFTs');
   end

   fp(2) = fopen(fpName2,'w+');
   if(fp(2) == -1)
      error('Could not open second temp file to hold FFTs');
   end

   fp(3) = fopen(fpName3,'w+');
   if(fp(3) == -1)
      error('Could not open third temp file to hold FFTs');
   end
end

c = normalizePow2(c);
CRamp = 1:length(c);
ExpM = degree : -1 : 0; % Modulator exponents
cPadded = zeros(1,padfactor*2^nextpow2(length(c)));

r = radii(Nradii); % Start with the largest search radius and work down.
cPadded(CRamp) = c .* r.^ExpM; % Fill in the first part; the rest is 0.0.

F2 = abs(rdft(cPadded,1));
nFreq = length(F2);
if(r > 1); % I found that this helps.
  F2 = F2 ./ (r^degree);
end
if(degree >= DEGREE_FOR_USING_FILES) 
   count = fwrite(fp(1), F2, 'double');
   clear F2;
end

r = radii(Nradii - 1);
cPadded(CRamp) = c .* r.^ExpM;
cPadded(CRamp) = c .* radii(Nradii - 1).^ExpM;
F3 = abs(rdft(cPadded,1));
clear cPadded
if(degree >= DEGREE_FOR_USING_FILES) 
   count = fwrite(fp(2), F3, 'double');
   clear F3;
end

LcPadded = padfactor*2^nextpow2(length(c));
SectorSize = pi/(nFreq-1); % To convert index to radians
j = Nradii - 2;
nFound = 0; % The total of roots found.
while(j)
  % fprintf('In grid search j=%d\n',j);
  % Cycle out F1; compute and cycle in the new F3.
  if(degree < DEGREE_FOR_USING_FILES) % If you are writing to files, F1 & F2 have been cleared.
    F1 = F2; 
    F2 = F3;
  end
  
  r = radii(j);
  cPadded = zeros(1, LcPadded);
  cPadded(CRamp) = c .* r.^ExpM;
  F3 = abs(rdft(cPadded,1));
  clear cPadded
  if(degree >= DEGREE_FOR_USING_FILES)
     count = fwrite(fp(3), F3, 'double');
     clear F3;
  end

  clear newIx; % Indices of new roots.
  if(degree < DEGREE_FOR_USING_FILES) 
    newIx = ix2dmin(degree, F1(:), F2(:), F3(:));
  else
    % Note, F1, F2, F3 files are rewound in ix2dminf()
    newIx = ix2dminf(degree, nFreq, fp(1), fp(2), fp(3));
  end
  
  % We have 2 goals: 1) Write new data to the oldest file. 2) Write to the third file in fp().
  % fp(1) is currently the oldest file so we want to write to it on the next pass. 
  % Thefore we cyclically permute fp() backwards so that the first file becomes the third file.
  if(degree >= DEGREE_FOR_USING_FILES)
    temp = fp(3);
    fp(3) = fp(1);
    fp(1) = fp(2);
    fp(2) = temp;
  end
  
  Newlength = length(newIx);
  if(Newlength)
    clear Theta Radii newRoots LnewRoots;
    Theta = SectorSize*(newIx - 1);
    Radii = radii(j+1) + zeros(Newlength,1);
    newRoots = Radii.*(cos(Theta) + i*sin(Theta));
    LnewRoots = length(newRoots);
    if(invert)
       rts1( (nRts1+1) : (nRts1+LnewRoots) ) = conj(1 ./ newRoots);
    else
       rts1( (nRts1+1) : (nRts1+LnewRoots) ) = newRoots;
    end     
    nRts1 = nRts1 + LnewRoots;
    nFound = nFound + LnewRoots;  
  end
  j = j - 1;
  
  if(DIAGNOSTICS == 2)
    % BATCH MODE: disabled diagnostic output
    % if(invert)
    %   fprintf('On iteration %d, radius=%14.9g the number found = %d\n', j, 1/radii(j+1), length(newRoots));
    % else
    %   fprintf('On iteration %d, radius=%14.9g the number found = %d\n', j, radii(j+1), length(newRoots));
    % end
  end
end

if(DIAGNOSTICS)
  % BATCH MODE: disabled diagnostic output
  % if(invert)
  %   fprintf('Grid searching radial range = %14.9g to %14.9g found %d putative roots\n', ...
  %           flip(1 ./ radii([1, Nradii])),nFound);
  % else
  %   fprintf('Grid searching radial range = %14.9g to %14.9g found %d putative roots\n', ...
  %           radii([1, Nradii]),nFound);
  % end
end

if(degree >= DEGREE_FOR_USING_FILES)
  fclose(fp(1)); fclose(fp(2)); fclose(fp(3)); pause(0);
  delete(fpName1); delete(fpName2); delete(fpName3);
else
  clear F1 F2 F3;
end

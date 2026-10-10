% function [radialM, thetaM] = findMultipliers(degree, radialMult, thetaMult);
%
% PURPOSE
% -------
% Determine thetaMult and radialM
%
% REQUIRED INPUT ARGUMENTS
% ------------------------
% degree - the polynomial degree
% radialMult - 1 or the value set on the command line
% thetaMult - 1 or the value set on the command line
%
% REQUIRED OUTPUT ARGUMENTS
% -------------------------
% radialM - The value to modify the default radial grid.
% thetaM - The multiplier applied to the polynomial length to determine FFT length.

% USES
% ----
% Nothing


% A 128-degree poly (i.e., 129 coefficients) uses twice as long a FFT as
% degree 127 (i.e. 128 coefficients). Thus the grid cells of a 128-degree
% are half the width in the theta direction of 127-degree polynomials. 
% Therefore, you get an AUMOMATIC boost from the longer FFT when degree
% reaches a power of 2. Therefore, any additional decreases using
% thetaMult should occur at midpoints between powers of 2, i.e. 1.5*2^n

% However, test at these midpoints. Perhaps you should instead increase 
% thetaMult at some other point.

% Yes indeed! degree=255 & thetaMult=1 requires a second pass 30% of
% the time on random coefficient smoothed by a Butterworth. You need to set 
% thetaMult=2 for degrees from 216 to 255 and then lower it back to 1
% after 255. After 255, you start using FFT's that are twice as long and 
% that is the equivalent of thetaMult=2 for degree < 255.

% If you are testing changes to deflateFFT() and you want to generate cases 
% that requires a LOT of deflated roots, then at degree==255 use thetaMult=1.
% The average is 10 roots by deflation, but I saw one polynomial that 
% required 9 deflations for 67 roots! This shows how powerful deflateFFT()
% is at doing multiple tries.

% Here are some runs that established the following values of thetaMult. I would like
% for the # of roots obtained by deflation to be < 10.

% degree	thetaMult	# of roots obtained	average
%							by deflation			seconds
% ------ ---------   -------------------	-------
%  28671       8     7,9,12                  138
%  28672      16     1,2,1                   146
%  32767      16     5,2                     172
%  32768       8     5,4                     287
%  36044       8     4,5,1                   299
%  65535      16     5,4                     649
%  72089       8     4,4,9                 1,218
% 114687      16     5,2,2                 2,267
% 131071      32     3,3,6                 2,733		
% 199,607     16     9
% 150,000     32     <= 5 roots on 5 tries
% 250,000     32     8,10,6,6              9,851
% 393,215     16     22                   29,540
% 393,217     32     11                   34,478
% 500,000     32     22,20                46,071,45234

% AUTHOR
% ------
% James W. Fox 12/28/2005

function [radialM, thetaM] = findMultipliers(degree, radialMult, thetaMult);

% Declare globals
global RAM DIAGNOSTICS MAX_ARRAY_LENGTH_EXPONENT;

radialM = 1;

if(degree < 50)
   thetaM = 2;
elseif(degree < 64)     % 2^6
  thetaM = 4;
elseif(degree < 70)     % 1.1*64
  thetaM = 2;
elseif(degree < 128)    % 2^7
  thetaM = 4;
elseif(degree < 141)    % 1.1*2^7
  thetaM = 2;
elseif(degree < 256)    % 2^8
  thetaM = 4;
elseif(degree < 282)    % 1.1*2^8
  thetaM = 2;
elseif(degree < 512)    % 2^9
  thetaM = 4;
elseif(degree < 563)    % 1.1*2^9
  thetaM = 2;
elseif(degree < 1024)   % 2^10
  thetaM = 4;
elseif(degree < 1126)   % 1.1*2^10
  thetaM = 2;
elseif(degree < 2048)   % 2^11
  thetaM = 8;
elseif(degree < 2253)   % 1.1*2^11
  thetaM = 4;
elseif(degree < 4096)   % 2^12
  thetaM = 8;
elseif(degree < 4505)   % 1.1*2^12
  thetaM = 4;
elseif(degree < 8192)   % 2^13
  thetaM = 8;
elseif(degree < 9011)   % 1.1*8192
  thetaM = 4;
elseif(degree < 16384)  % 2^14
  thetaM = 8;
elseif(degree < 18022)  % 1.1*2^14
  thetaM = 4;
elseif(degree < 28672)  % 1.75*2^14
  thetaM = 8;
elseif(degree < 32768)  % 2^15
  thetaM = 16;
elseif(degree < 36045)  % 1.1*2^15
   thetaM = 8;
elseif(degree < 65536)  % 2^16
  thetaM = 16;
elseif(degree < 72090)  % 1.1*2^16
  thetaM = 8;
elseif(degree < 114688) % 1.75*2^16
   thetaM = 16;
elseif(degree < 131072) % 2^17       Thoroughly checked through here
  thetaM = 32;
elseif(degree < 144179) % 1.1*2^17
  thetaM = 16;
elseif(degree < 262144) % 2^18
   thetaM = 32;
elseif(degree < 288358) % 1.1*2^18 Two tries at 288356 & 16 needed 4 & 7 roots by deflation.
   thetaM = 16;
elseif(degree < 350000) % 400,000 & 32 needed 16 by deflation. 350,000 & 32 needed 8.
   thetaM = 32;       
elseif(degree < 524288) % 2^19
  % thetaM=64 gives "Out of memory" for 2 GB RAM on Windows but works for 2 GB. on Linux.
  % 500,000 & thetaM=32 & radialM=1 needed 40 roots by deflation. Dangerous!
  % 500,000 & thetaM=32 & radialM=2 needed 5 roots by deflation.
  % 500,000 & thetaM=32 & radialM=4 needed 1 root by deflation.
  if(RAM <= 2) % Works well with RAM=1 & RAM=2 - missed 1 upper half-plane root.
    thetaM = 32; % Even with the recent changes, 64 ran out of memory after 2 inner loop iterations in grid search.
    radialM = 2;
  else              
    thetaM = 64;
  end
elseif(degree < 1048576) % 2^20
  if(RAM <= 1)
    thetaM = 16;
    radialM = 8; % Tested once. Grid search missed 11 UHP roots. Oddly, radialM=16 also missed 11.
  elseif(RAM <= 2)
    thetaM = 32;
    radialM = 2; % "1" failed once and worked poorly 2 other times. "2" worked great 3 times out of 3.
  else
    thetaM = 64;
    radialM = 1; % Tested once. Grid search missed 1 upper half-plane root. Thanks Ryan King!
  end
elseif(degree < 2097152) % 2^21 
  if(RAM <= 2)
    thetaM = 16;  % Tried on degree 2,080,000. Using radialM=4 missed 21 upper half-plane roots.
    radialM = 32; % RadialM=16 missed 17. RadialM=32 missed 10. All worked. Thanks Sid Burrus!
  elseif(RAM <=4)
    thetaM = 32;  % Tried on degree 2e6 under Linux.
    radialM = 4;  % Grid search found every last root! Thanks Ryan King!
  elseif(RAM <= 8)
    thetaM = 64; % This has never been tested.
    radialM = 2;
  else
    thetaM = 128; % This has never been tested.
    radialM = 1;
  end
else  % This has not been tested because it takes too long - estimated 16 cpu days.
  if(RAM <= 2)
    error('You do not have enough RAM');
  elseif(RAM <= 4)
    thetaM = 16;
    radialM = 32;
  elseif(RAM <= 8)
    thetaM = 32;
    radialM = 8;
  else
    thetaM = 64;
    radialM = 2;
  end   
end             

radialM = radialM * radialMult;
thetaM = 2^nextpow2(thetaM * thetaMult); % Apply thetaMult and make sure it is a power of 2.

if(thetaM*2^nextpow2(degree+1) >  2^MAX_ARRAY_LENGTH_EXPONENT)
  error('You do not have enough memory to factor this high a degree.');
end

if(DIAGNOSTICS)
  % BATCH MODE: disabled diagnostic output
  % fprintf('\nthetaM=%d \tradialM=%d \tRAM=%d        %s, lroots version 2.\n', thetaM, radialM, RAM, datestr(now));
  % hwrite(['thetaM=',num2str(thetaM),', radialM=',num2str(radialM),', RAM=',num2str(RAM),'   ', datestr(now),', lroots version 2.']);
end

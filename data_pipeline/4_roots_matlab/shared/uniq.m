% function [uniqueZ, zErr, dupIndicies] = uniq(z, zErr);
%
% PURPOSE:
%   Find the unique values of a fuzzy complex array. If two values differ by at
%   most the error estimate, declare them equal and eliminate the duplicate.
%
% REQUIRED INPUT:
%   z -  Any 1D complex array.
%   zErr - An array or error estimates for z. If only a single number is supplied,
%          it is repeated. If [] is supplied 1e-7 is used.
 %
% REQUIRED OUTPUT:
%   uniqueZ - The unique values of z.
%
% OPTIONAL OUTPUT
%   zErr - An array or error estimates for z. If [] is supplied, an array of 1e-7's
%          is used. If only a single number is supplied, it is repeated. 
%   dupIndicies - Indicies of duplicate elements. 
%                 length(z) = length(uniqueZ) + length(dupIndicies)

% 	WARNING
%  The error extimates may be too small. I have found it necessary to multiply
%  the errorestimates by 10 before calling uniq. Cf lroots and ????

% USES:
%   nothing
%
% AUTHOR:
%   James W. Fox 9/24/2001. 

% CHANGES
%   Gary Sitton provided an excellent suggestion for speeding up the original version.

function [uniqueZ, zErr, dupIndicies] = uniq(z, zErr);

DIAGNOSTICS = 1; % 0=none, 1=minimal, 2=full diagnostic messages

if(length(z) <= 1) % User supplied only 0 or 1 value.
  uniqueZ = z;
  dupIndicies = [];
  return
elseif(nargin ~= 2)
   error('Supply exactly 2 arguments.');
elseif(isempty(zErr))
   zErr = 1e-7 + zeros(size(z));
elseif(length(zErr) == 1)
   zErr = abs(zErr) + zeros(size(z));
else                          
  if(length(zErr) ~= length(z))
    error('length($1) must equal length($2).')
  elseif (min(zErr) < 0)
    error('$2 must be not contain negative values.')
  end
end

% Sort by complex argument
ang = angle(z);
[angZ, iangZ] = sort(ang);
z = z(iangZ);            % angZ and iangZ are also used below.
zErr = zErr(iangZ);

Lz = length(z);
j = 1;
while(j < Lz)
%   if(DIAGNOSTICS)
%     if(~mod(j,1000))
%        fprintf(1,'\t uniq() on root # %d\n', j);
%     end
%   end
   
   % At the very end of the program I collect all the unique values.
   % Alternatively, save z(j) now; it will be a unique value.
   
   % Consider z(j) as the start of a cluster of roots with arguments
   % between ang(z(j)) and ang(z(j))+zErr(j). (Usually the "cluster"
   % will only contain 1 value.) Scan the cluster for duplicates
   % of z(j) and flag them by resetting them to NaN.
   
   % It is possible to have another root in the cluster with 
   % essentially identical magnitude but very different argument.
   % Take note if there are any such false duplicates in the cluster.
   
   % The original code tested the difference to see if it was < zErr(j)
   % which defaulted to 1e-7. I since switched to using the final corrections
   % as an estimate of the error. This caused a problem.  Here is an example
   % from a 1000-degree polynomial when uniq failed and here is why:

   %» rts(52:53)
   % ans = -4.650444723145865e-001 +8.802772363578753e-001i
   %       -4.650444723145864e-001 +8.802772363578753e-001i
   %
   %» rtsErr(52:53)
   % ans = 1.044048714879764e-051  <======== Wildly too small.
   %       8.462308737430829e-015
   %
   % The difference between these two roots was 5e-17 but that was not less than
   % rtsErr(52)=1e-51 so they were declared to be different roots. BTW, I 
   % multiplied the actual error estimates by 200 in lroots.
   
   % I am changing the test to zErr(j)+zErr(k). This works in the above situation
   % because |z(j)-z(k)| = 5e-17 < zErr(j)+zErr(k) = 1e-51 + 8e-15.
   
   % BTW, if you have an ill-conditioned root and zErr(k) is huge, then uniq may
   % unnecessarily eliminate a lot of roots. I think what you want it is:
   % z(j) and z(k) are false double roots if |z(j)-z(k)|<min(zErr(j),zErr(k)) but
   % both error estimates must be accurate. You can not have 1e-51 as an error 
   % estimate.
   
   falseDuplicate = 0; % Are there any?

   k = j + 1; % Index used to scan the rest of the cluster.
   while(k <= Lz & (angZ(k)-angZ(j)) <= zErr(j)+zErr(k)) % If z(k) is in the cluster:
      if(~isnan(z(k)))
         if(abs(z(k)-z(j))) <= zErr(j)+zErr(k)
            z(k) = NaN;
         else
            falseDuplicate = 1; % z(k) is truly different from z(j).
         end
      end
      k = k+1;
   end
   
   % Where shall we start the next cluster? At this point k is the index
   % of the first value after the end of the last cluster. However, if 
   % the last cluster contained any value that was truly different from 
   % z(j), then start the new cluster at that value.
   
   if(falseDuplicate) % Start at the different value inside the cluster.
      j = j + 1;
      while(j<Lz & isnan(z(j))) % Skip past Nan's
         j = j + 1;
      end
   else               % Start at the next point after the last cluster.
      j = k;
   end
end

k = find(~isnan(z)); % Indicies of unique values.
uniqueZ  = z(k);

if(nargout > 1)
   zErr = zErr(k);
end

if(nargout > 2)
   dupIndicies = iangZ(find(isnan(z)));
end

% Diagnostic message:
% fprintf('\n\t\t\t\t\t\t<<uniq is throwing away %d values.>>\n',length(dupIndicies));

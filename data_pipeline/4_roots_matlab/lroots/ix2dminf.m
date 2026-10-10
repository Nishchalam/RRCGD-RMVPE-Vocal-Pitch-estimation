% function ixRelMin = ix2dminf(degree, nFreq, fp1, fp2, fp3)

% PURPOSE
%   Search for 2D relative minima. Cf ix2dmin. This assumes the data has been written
%   to 3 files with file pointers fp1, fp2, fp3.
%
% REQUIRED INPUT:
%   degree - The degree.%   fp1 - file pointer to the first file containing data.
%   nFreq - length of the arrays in the files.
%   fp2 - file pointer to the second file containing data.
%   fp3 - file pointer to the third file containing data.
%
% OUTPUT:
%   ixRelMin - indices of relative 2D minima. A root is likely to be nearby.

% USES:
%   nothing
%
% AUTHOR:
%   James W. Fox 12/19/05

function ixRelMin = ix2dminf(degree, nFreq, fp1, fp2, fp3)

bucketSize = 2^26; % How many values to read at once. Use 2^26 for 1/2 GB. Use 2^6 for tests with degree 1000.

bucketSize = round(nFreq/4); % ????

bucketSize = round(min(bucketSize, nFreq));
nBuckets = floor(nFreq/bucketSize); 

frewind(fp1); frewind(fp2); frewind(fp3);
ixRelMin = zeros(1,degree);
nSaved = 0; % Number currently saved in ixRelMin.

for m = 1:nBuckets
  if(m < nBuckets)   
    F1 = fread(fp1,bucketSize,'double');
    F2 = fread(fp2,bucketSize,'double');
    F3 = fread(fp3,bucketSize,'double');
    fseek(fp1, -16, 0); % Rewind 2 samples so that on the next read we will get a full 
    fseek(fp2, -16, 0); % complement of points around the last sample of the last read.
    fseek(fp3, -16, 0);
  else % Read to EOF,
    F1 = fread(fp1,inf,'double');
    F2 = fread(fp2,inf,'double');
    F3 = fread(fp3,inf,'double');
  end
  
  ixBucketBegin = bucketSize*(m-1) - 2*m + 3; % Indices in the full array.
  ixBucketEnd   = m*(bucketSize - 2) + 2; % Add bucketSize-1 to the previous line.
  ix = ix2dmin(degree, F1, F2, F3) + ixBucketBegin - 1;  % Indices in the full array.
  clear F1 F2 F3;
  
  % If ixBucketBegin or ixBucketEnd is a relative minimum
  % (and not the first or last sample of the entire array) then eliminate it.
  if(length(ix))
    if(ix(1) == ixBucketBegin & m > 1)
      ix = ix(2:end);
    elseif(ix(1) == ixBucketEnd & m < nBuckets)
      ix = [];
    end
    if(length(ix) > 0 & ix(end) == ixBucketEnd & m < nBuckets)
      ix = ix(1:end-1);
    end
    Lix = length(ix);
    ixRelMin(nSaved + 1: nSaved + Lix) = ix;
    nSaved = nSaved + Lix;
  end
  clear ixBucketBegin ixBucketEnd ix Lix
end

frewind(fp1); frewind(fp2); frewind(fp3); % So they can be written to at the beginning.
ixRelMin = ixRelMin(1:nSaved)';

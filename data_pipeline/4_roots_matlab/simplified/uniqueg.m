% Find an array's unique values. If magnitudes of
% two values differ by at most gap, eliminate the 
% duplicates; all values = cluster's first value.
%
% Input:
%   z = Any complex vector.
%
% Output:
%   uz = Clustered values of z.
%
% Uses:
%   sort.
%
% uz = uniqueg(z);
%
function uz = uniqueg(z);
%
[absZ, iabsZ] = sort(abs(z)); z = z(iabsZ); % Sort.
Lz = length(z); uz = zeros(Lz, 1); % Allocate space.
gap = 1e-7; uzi = 0; j = 1; % Initialize.
while (j <= Lz) % Scan for duplicates.
    k = j + 1; % Index inside cluster.
    falseDup = 0; % Clear flag.
    while ((k <= Lz) & ...
            (absZ(k) - absZ(j)) <= gap) % Scan the cluster.
        if (~isnan(z(k))) % Skip deleted values.
            if (abs(z(k) - z(j)) <= gap) % Cluster member?
                z(k) = NaN; % Delete value.
            else
                falseDup = 1; % Truly different.
            end
        end
        k = k + 1; % Next cluster index.
    end
    uzi = uzi + 1; uz(uzi) = z(j); % Save next unique value.
    if (falseDup) % Mixed cluster?
        j = j + 1;
        while ((j < Lz) & isnan(z(j))) % Skip deleted values.
            j = j + 1;
        end
    else
        j = k; % Start after the cluster.
    end
end
uz = uz(1:uzi); % Truncate array.



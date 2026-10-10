% Ping-pong root shuffle/sort.
%
% Input:
%   r = column vector of angle sorted roots.
%
% Output:
%   sr = ping-pong sorted column vector.
%
% sr = pprsort(deg, rmin);
%
function sr = pprsort(r);
%
rn = length(r); r2 = 2*fix(rn/2); % Get lengths. 
sr = [r(1:2:r2-1).'; r(r2:-2:2).']; % Ping-pong rows.
sr = sr(:); % Merge to rows.
% Last element?
if (length(sr) < rn)
    sr(rn) = r(rn);
end
% Done.

% Radial distribution generator.  Generate
% radial distribution values from 1 to rmin.
%
% Input:
%   deg = degree of (random) polynomial,
%   rmin = minimum radius value < 1.
%
% Output:
%   rd = ordered radial distribution values.
%
% rd = radiald(deg, rmin);
%
function rd = radiald(deg, rmin);
%
inc = 1/8; rad = 1; rn = 0; % Initialize.
rd = zeros(deg, 1); % Allocate space.
while (rad >= rmin | rn < 3) % Fill radius array.
    rn = rn + 1; rd(rn) = rad; % Store next value.
    if (rad >= (1 - 16/deg)) % Test grid radius.
        rdt = 0.012*deg*(rad - 1 - 8.09/deg)^2; % Near to u.c.
    else
        rdt = min(1, (1 - rad - 10/deg)); % Far from u.c.
    end
    rad = rad - inc*rdt; % Next radius.
end
rd = rd(1:rn); % Truncate array.
% Done.
% Use the FFT for a fast grid search to find 
% approximate roots; use newton to polish them.
%
% Input:
%   c = Polynomial coefficients in decending order.
%   radii = Vector of search radii.
%   padfact = Increasing this value makes the search 
%           grid finer in the theta direction.
%   invert = If ~= 0, flip the coeffs. & invert roots. 
%            
% Output:
%   root = A vector of polished roots found by searching
%           the radial range.
% Uses:
%   newton, uniqueg, fft.
%
% root = search(c, radii, padfactor, invert);
%
function root = search(c, radii, padfact, invert);
%
c = c(:); cn = length(c); Nradii = length(radii); % Lengths.
Degree = cn - 1; root = zeros(Degree, 1); % Allocate space.
ExpM = [Degree:-1:0]'; % Modulator exponents.
% Flip coefficients?
if (invert) 
    c = c(cn:-1:1); 
end
n = 2^ceil(log2(padfact*cn)); n2 = n/2 + 1; % Next power-of-2.
k = [2:n2+1]; kp1 = k + 1; km1 = k - 1; % Index arrays.
FA = abs(fft(c.*(radii(Nradii).^ExpM), n)); % F2 radius. 
F2 = [FA(2); FA(1:n2); FA(n2-1)]; % Pad F2.
FA = abs(fft(c.*(radii(Nradii-1).^ExpM), n)); % F3 radius.
F3 = [FA(2); FA(1:n2); FA(n2-1)]; % Pad F3.
SectorSize = pi/(n2 - 1); % Index in radians.
j = Nradii - 2; nr = 0; % Initialize.
while (j) % Scan radial grid.
    F1 = F2; F2 = F3; % New F1 & F2.
    FA = abs(fft(c.*(radii(j).^ExpM), n)); % New F3.
    F3 = [FA(2); FA(1:n2); FA(n2-1)]; F2k = F2(k); % Middle points.
    NewIx = find((F2k < F1(km1)) & (F2k < F3(km1)) & (F2k < F1(k)) ...
        & (F2k < F3(k)) & (F2k < F1(kp1)) & (F2k < F3(kp1)) ... 
        & (F2k < F2(km1)) & (F2k < F2(kp1))); % Vectorized minima test.
    lx = length(NewIx);
    if (lx) % Any minima found?
        Theta = SectorSize*(NewIx - 1); % Get root angles.
        root(nr+1:nr+lx) = radii(j+1)*complex(cos(Theta), sin(Theta));
        nr = nr + lx; % Increment count.
    end
    j = j - 1; % Increment radial index.
end
if (nr)  % Any roots found? 
    root = uniqueg(newton(c, root(1:nr), 0)); % Polish roots.
    % Invert the roots?
    if(length(root) & invert) 
        root = conj(1./root); 
    end
else
    root = []; % No roots.
end
% Done.
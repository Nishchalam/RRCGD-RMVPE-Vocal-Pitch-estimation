% Use roots of a real polynomial to 
% unwrap its Fourier phase spectrum.
% 
% Input:
%   x = Real polynomial coefficients, 
%       highest order term 1st.
%   ep = Optional leading zero relative magnitude 
%       threshold, default = 1e-12.
%   z = Optional complex root input.
%
% Ouput:
%   ph = Unwrapped positive phase spectrum in radians,
%   err(1) = Phase error from unwrapper,
%   err(2) = Coeff. error from root finder,
%   amp = Unwrapped positive amplitude spectrum,
%   z = Half-plane complex polynomial roots.
%
% [ph, err, amp, z] = runwrap(x, ep, z);
%
% Uses: simp_fl.
%
function [ph, err, amp, z] = runwrap(x, ep, z)
% Default ep?
if (nargin < 2)
    ep = 1e-12;
end
% Normalize coeffs.
x = x(:)'/max(abs(x));
% Find half-plane roots?
if (nargin < 3 | isempty(z))
    [z, er, in] = simp_fl(x, ep); err(2) = er;
else
    in = find(abs(x) > ep);
    err(2) = 0;
end
% Convert to half-plane.
z = z(find(imag(z) >= 0));
% Truncate leading zeros.
N = length(x); x = x(in(1):N);
% Initialize counters, etc.
NCxMin = 0; NUCR = 0; N_1 = 0; N2 = fix(N/2)+1; 
ph = zeros(1, N2); amp = abs(fft(x));
% DFT cos/sin table.
k = [0:N2-1]*(2*pi/N); cs = complex(cos(k), -sin(k));
% Main loop
for j = 1:length(z)
    % Get next root.
    z1 = z(j); r = abs(z1); ir = imag(z1);
    % Very near the unit circle?
    if (abs(1 - r) < ep)
        if (ir == 0) % Real unit cirle root.
            if (abs(1 - z1) < ep) 
                N_1 = N_1 + 1; % z1 == 1
            end
        else % Complex unit circle root.
            NUCR = NUCR + 1;
            ix = [(1 + ceil(N*angle(z1)/(2*pi))):N2];
            amp(ix) = -amp(ix); % Unwind amplitude spectrum.
        end
        % Unit circle counter.
        NUCR = NUCR + 1;
    else
        % Outside the u.c.?
        if (r > 1)
            z1 = 1/z1; ir = imag(z1);
        end
        % Root components.
        rr = real(z1); 
        % Compute DFT of root.
        if (ir == 0) % Real case.
            spec = 1 - cs*rr;
        else % Complex case.
            spec = (1 - cs.*(2*rr - cs*(rr*rr + ir*ir)));
        end
        % Add/subtract phase.
        if (r < 1) % Inside the u.c.
            ph = ph + angle(spec);
        else % Outside the u.c.
            ph = ph - angle(spec);
            % Min. Phas. counter.
            if (ir == 0) % Real root.
                NCxMin = NCxMin + 1;
            else % Complex root.
                NCxMin = NCxMin + 2;
            end
        end
    end
end
% Add linear corrections.
ph = ph - pi*((2*NCxMin + NUCR)*[0:1/N:.5] - N_1*.5);
% D.C. phase.
if (amp(1) < ep)
    ph(1) = 0; 
end
% Reconstruct the signal.
amp = amp(1:N2); y = amp.*complex(cos(ph), sin(ph));
% Force conjugate symmetry.
y = real(ifft([y conj(y(N2-rem(N+1, 2):-1:2))]));
% Scale correction and error.
nm = (y*x')/(y*y'); amp = nm*amp; err = max(abs(x - nm*y));
% Done.


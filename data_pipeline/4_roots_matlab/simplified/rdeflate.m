% Divide a real polynomial c by a set of half-plane 
% complex roots r to form the deflated real polynomial 
% coefficients rc.  Uses frequency domain via FFT/DFT/
% IFFT to obtain the deflated coefficient set.
%
% Input:
%    c = coefficient array,
%    r = complex root array,
%
% Output:
%    rc = reduced coefficient array,
%    r1 = 1st quotient remainder term.
%
% Uses: 
%   fft, ifft.
%
% [rc, r1] = rdeflate(c, r);
%
function [rc, r1] = rdeflate(c, r);
%
cn = length(c); rn = length(r); % Lengths.
m = 2^(ceil(log2(cn))); % FFT length.
m2 = m/2 + 1; % Nyquist index.
fc = fft(c, m); % FFT of coeffs.
fc = fc(1:m2); % Half-plane only.
fc = fc(:); % Force to column.
k = [0:m2-1]'*(2*pi/m); % DFT cos/sin table.
cs = complex(cos(k), -sin(k));
i = 1; rrn = 0; % Initialize.
while (i <= rn) % Division loop.
    rr = real(r(i)); ir = imag(r(i)); % Components.
    if (ir ~= 0) % Complex root?
        fc = fc./(1 - cs.*(2*rr - cs*(rr*rr + ir*ir)));
        rrn = rrn + 1;
    else % Real root.
        fc = fc./(1 - cs*rr);
    end
    i = i + 1; rrn = rrn + 1; % Increment.
end
rc = ifft([fc; conj(fc(m2-1:-1:2))]); % Invert.
r1 = rc(cn-rrn+1); % Remainder.
rc = real(rc(1:cn-rrn)); % Truncate.

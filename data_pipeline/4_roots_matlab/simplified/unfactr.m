% Unfactor a set of complex half or full-plane  
% roots into a set of real polynomial coefficients 
% uc.  Uses the frequency domain via DFT/IFFT to
% obtain the unfactored coefficient set.
%
% Input:
%   r = half or full-plane complex root array.
%
% Output:
%   uc = unfactored coefficient array.
%
% Uses:
%   ifft.
%
% uc = unfactr(r);
%
function uc = unfactr(r);
% 
rn = length(r); % Array length.
m = 2^(ceil(log2(2*rn+2))); % FFT length.
m2 = m/2 + 1; % Nyquist index.
k = [0:m2-1]'*(2*pi/m); % DFT cos/sin table.
cs = complex(cos(k), -sin(k));
i = 1; rrn = 0; fc = ones(m2, 1); % Initialize
while (i <= rn) % Multiplication loop.
    rr = real(r(i)); ir = imag(r(i)); % Components.
    if (ir ~= 0) % Complex root?
        fc = fc.*(1 - cs.*(2*rr - cs*(rr*rr + ir*ir)));
        rrn = rrn + 1;
        % Conjugate pair?
        if (i ~= rn & ir == -imag(r(i+1)))
            i = i + 1;
        end
    else % Real root.
        fc = fc.*(1 - cs*rr);
    end
    i = i + 1; rrn = rrn + 1; % Increment.
end
fc = ifft([fc; conj(fc(m2-1:-1:2))]); % Invert
uc = real(fc(1:rrn+1)); % Truncate.
% Done.
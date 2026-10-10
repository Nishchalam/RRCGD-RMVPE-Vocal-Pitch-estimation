% Newton's method optimized for large real polynomials. 
% If no converge in maxiter steps, no root is returned.  
% Coefficients must be properly normalized and truncated. 
%
% Input:
%   x = desending order real polynomial coefficents,
%   zin = an array of complex half-plane starting roots,
%   fpflg = full-plane roots flag (default = 1), 
%           1 => full-plane roots supplied
%           0 => half-plane roots supplied.
%   maxit = Maximum number of iterations,
%   mindz = Minimum required correction,
%
% Output:
%   zout = Vector of polished half-plane roots.
%
% Uses:
%   polyval.
%
% zout = newton(x, zin, fpflg, maxit, mindz);
%
function zout = newton(x, zin, fpflg, maxit, mindz);
% Default fpflg?
if (nargin < 3) 
    fpflg = 1; 
end
% Convert to half-plane input?
if (fpflg) 
    zin = zin(find(imag(zin) >= 0)); 
end
x = x(:); nx = length(x); nzin = length(zin); % Lengths.
deg = nx - 1; zout = zeros(nzin, 1); % Allocate spcae.
% Default maxit? 
if (nargin < 4) 
    maxit = 32; 
end
% Default mindz?
if (nargin < 5) 
    mindz = 2^(-46); 
end
k = 0; dx = []; flipped_x  = []; dflipped_x = []; % Initialize.
for j = 1:nzin % Polish all roots.
    z = zin(j); ntries = 0; finish = 0; flipflg = -1; % Initialize.
    while (finish < 2) % Polish current root.
        ntries = ntries + 1; % Increment tries.
        % Divergence?
        if (ntries > 1 & isinf(z)) 
            finish = 2; absdz = 1.0; % Reject root.
        else
            % Outside u.c.? Invert root/flag. 
            if (abs(z) > 1) 
                z = 1/z; flipflg = -flipflg; 
            end
            if (flipflg == -1) % Normal Newton?
                if (isempty(dx)) % First time?
                    dx = [deg:-1:1]'.*x(1:deg); % Derivative.
                end
                f  = polyval(x, z); % Polynomial value.
                df = polyval(dx, z); % Derivative value.
            else % Inverted Newton.
                if (isempty(flipped_x)) % First time?
                    flipped_x  = x(nx:-1:1); % Flipped versions.
                    dflipped_x = [deg:-1:1]'.*flipped_x(1:deg);
                end
                f  = polyval(flipped_x, z); % Flipped value.
                df = polyval(dflipped_x, z); % Flipped derivative.
            end
            % Valid correction?
            if (f == 0) 
                dz = 0; % Zero update.
            else 
                dz = f/df; % Newton update.
            end 
            absdz = abs(dz); z = z - dz; % Update root.
            % Convergence?
            if (absdz < 1e-16) 
                finish = 2; % Really finished.
            elseif ((absdz < mindz) | (ntries >= maxit))
                finish = finish + 1; % Keep going.
            end  
            % Force to real?
            if ((ntries > 4) & (abs(imag(z)) < 1E-13)) 
                z = real(z); end
        end
    end % End iteration for this root.
    % Force to real?
    if (abs(imag(z)) < 1E-7) 
        z = real(z); 
    end  
    % Accept the root?
    if (((absdz < max(1E-7, mindz)) | ((absdz < 1e-3) & ...
       (abs(f) < 1e-13) & abs(z) > 0))) 
        % Invert the root?
        if (flipflg == 1) 
            z = 1/z; 
        end
        k = k + 1; % Force upper half-plane.
        zout(k) = complex(real(z), abs(imag(z)));
    end
end % End main root iteration loop.
zout = zout(1:k); % Truncate root array.
% Done.
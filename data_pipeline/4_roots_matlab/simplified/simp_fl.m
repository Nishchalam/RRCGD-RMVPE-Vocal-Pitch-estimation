% Simple L-F program.  A minimal, robust Matlab
% implementation of the Lindsey-Fox algorithm.
%
% Input:
%   c = descending order real polynomial 
%       coefficients, degree d = length(c)-1.
%   t = leading/trailing coefficient realative
%       threshold, default = 1e-12.
%   p = print diagnostic information if ~= 0,
%       default = 0.
%
% Output:
%   r = d full-plane complex polynomial roots,
%   e = maximum "unfactored" coefficient error.
%   in = array of index bounds of "nonzero" coeffs.
%       in(1) = first used coeff. index,
%       in(2) = last used coeff. index.
%
% Uses:
%   newton, uniqueg, radiald, unfactr, pprsort,
%   deflate, search, deconv, and roots.
%
% [r, e, in] = simp_fl(c, t, p);
%
function [r, e, in] = simp_fl(c, t, p);
% Default p?
if (nargin < 3) 
    p = 0; 
end
% Default t?
if (nargin < 2) 
    t = 1e-12; 
end
% Normalize coefficients.
c = c(:)/max(abs(c(:))); ocn = length(c);
% Any leading/trailing zeros?
i = find(abs(c) > t); in(1) = i(1);
i = find(abs(c(ocn:-1:in(1))) > t); 
% Truncate off zeros.
in(2) = ocn-i(1)+1; c = c(in(1):in(2)); cn = length(c);
if (p)
    disp(sprintf('Active index range: %d to %d.', in(1), in(2)));
end
% Bypass F-L algorithm?
if (cn < 64)
    % Polish half-plane roots.
    r = newton(c, roots(c));
else
    % Grid search parameters
    tm = max(1, sqrt(cn)/10); rmin = 1/(1 + 256/cn);
    if (p) 
        disp(sprintf('Radial limits: %e, %e', rmin, 1/rmin)); 
    end
    % 1 >= rd >= rmin.
    rd = radiald(cn-1, rmin);
    % Roots inside/on u.c.
    r = search(c, [1/rd(2); rd], tm, 0);
    if (p) 
        disp(sprintf('%d half-plane roots <= 1.', length(r))); 
    end
    % Roots outside unit circle.
    rr = search(c, rd, tm, 1);
    if (p) 
        disp(sprintf('%d half-plane roots > 1.', length(rr))); 
    end
    % Combine all roots.
    r = uniqueg([r; rr]);
    % Full-plane root count.
    nr = 2*length(r) - length(find(imag(r) == 0));
    if (p) 
        disp(sprintf('%d full-plane roots found.', nr)); 
    end
    % Ping-pong the roots.
    r = pprsort(r); 
    % Deflate polynomial?
    if (nr < cn-1)
        % Separate roots near u.c.
        uc = (abs(abs(r) - 1) < 1e-12);
        rr = r(find(uc)); cr = r(find(~uc));
        % Divide out normal roots.
        if (p) 
            disp('Deflating the polynomial...'); 
        end
        [rc, e] = rdeflate(c, cr);
        if (p) 
            disp(sprintf('%d roots left, error = %e', ...
                length(rc)-1, abs(e))); 
        end
        % Any roots near u.c.?
        if (~isempty(rr))
            % Divide out u.c. roots.
            [rc, e] = deconv(rc, unfactr(rr));
            if (p) 
                disp(sprintf('U.C. deflation error = %e', ...
                    abs(e(1)))); 
            end
        end
        % Polish half-plane roots.
        rr = newton(c, roots(rc), 1, 5);
        % Add residual roots.
        r = [r; rr];
    end
end
% Unfactor all the roots.
cc = unfactr(r);
if (p) 
    disp(sprintf('%d out of %d nonzero roots + %d roots at zero.', ...
        length(cc)-1, length(c)-1, ocn-in(2))); 
end
% Compute coefficient error?
if (length(c) >= length(cc))
    cc = [cc; zeros(length(c)-length(cc), 1)];
else    
    if (p) 
        disp(' Too many roots found.'); 
    end ; 
    e = Inf; r = []; % Error. 
    return;
end
% Normalized error.
cc = cc*((cc'*c)/(cc'*cc)); e = abs(max(c - cc));
if (p) 
    disp(sprintf('Maximum coefficient error = %e', e)); 
end
% Separate complex & real roots.
cr = (imag(r) ~= 0);
rr = r(find(cr)); r = r(find(~cr));
r = [rr; rr; r; zeros(ocn-in(2), 1)];
% Pair conjugate roots.
cn = 2*length(rr);
r(1:2:cn-1) = rr; r(2:2:cn) = conj(rr);
% Done.
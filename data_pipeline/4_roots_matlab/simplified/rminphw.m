% Find the minimum phase factoring of a wavelet
% using polynomial root finding and inversion.
% Coeffs. ordered in increasing powers of 1/z.
% 
% Input:
%   x = vector of real polynomial coeffs.
%
% Output:
%   mx = minimum phase equivalent of x.
%   zo = roots of original polynomial,
%   zm = roots of min. phas. polynomial.
%
% [mx, zo, zm] = rminphw(x);
%	
function [mx, zo, zm] = rminphw(x);
% Find all the roots.
[zo, e] = simp_fl(x);
% Error?
if (e > 1e-6)
    error(' Error in root finder!');
end
% Find the roots outside.
i = find(abs(zo) > 1);
% Invert these roots.
zm = zo; zm(i) = conj(1./zo(i));
% Build new polynomial.
mx = unfactr(zm);
% Rescale output signal. 
mx = mx*(norm(x)/norm(mx));
% Change the sign?
if mx(1)*x(1) < 0
    mx = -mx;
end
% Add any zeros?
mx = [mx; zeros(length(x)-length(mx), 1)];
% Done.


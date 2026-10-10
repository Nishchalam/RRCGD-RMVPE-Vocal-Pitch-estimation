% Plot roots of a polynomial 
% in the complex plane.
% p = real polynomial coeffs.,
%     or complex roots.
% ps = plot symbol string.
%
% proots(p, ps);
%
function proots(p, ps)
% Default symbol
if (nargin < 2)
   ps = 'r.';
end
% Get poly. size
[m, n] = size(p);
% Column vector?
if (n == 1)
    p = p';
    [m, n] = size(p);
end
% Complex roots?
if (any(imag(p(:))) ~= 0)
    rp = p(:);
else
    % Space for roots
    rp = zeros(n-1, m);
    % Loop for roots
    for i = 1:m
        rp(:, i) = simp_fl(p(i, :));
    end
end
% Plot unit circle
th = (2*pi/100)*[0:100]';
plot(cos(th), sin(th), 'b-');
axis('square');
hold on;
plot(0, 0, 'b+');
% Plot all roots
for i = 1:m
    % Test for linear section
    if (n == 3) & (p(i, 3) == 0)
        plot(real(rp(2, i)), imag(rp(2, i)), ps);
    else
        plot(real(rp(:, i)), imag(rp(:, i)), ps);
    end
    plot(0, 1.1, ' ');
    hold on;
end
% Draw axes
plot([-1.1 1.1], [0 0], '--');
plot([0 0], [-1.1 1.1], '--');
% Scale axes
axis('equal');
title('Complex Polynomial Roots');
xlabel('Real Axis');
ylabel('Imaginary Axis');
hold off;
%


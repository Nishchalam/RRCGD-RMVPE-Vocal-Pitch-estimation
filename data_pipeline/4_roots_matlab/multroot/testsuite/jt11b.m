function [p,z] = jt11b
%
%  test polynomial suggested by Jenkins and Traub
%
    m = 20;
    y = exp(i*pi*[(1-m):(m-1)]/(2*m));
    z = 0.9*exp(i*pi*[m:3*m]/(2*m));
    u = [y,z];
% JWF
%   p = poly(u);
    p = poly(leja(u));

    z = [[y';z'],ones(4*m,1)];
    fprintf('\n');
   z = displayRoots(z,2); % This also clips small parts and sorts
    %fprintf(' Roots are all simple and sensitive, see plot ');
    %plot(real(z(:,1)),imag(z(:,1)),'+');
    % JWF. No the roots are not sensitive. Instead, the coefficients are sensitive to the
    % order in which the roots are multiplied together. Doing a leja sort before
    % multiplying them together helps immensely.      
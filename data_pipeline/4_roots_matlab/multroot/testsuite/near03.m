function [p,z] = near03
%
%  test polynomial suggested by Z. Zeng
%  It is suggested that MultRoot input gamma = 10.
%
    e = 0.001;
% JWF
% Using leja helps immensely. It still fails sometimes, but works much more often.
%   p = poly([(1-e)*ones(1,20),ones(1,20),-0.5*[1,1,1,1,1]]);
%   z = [1-e, 20; 1, 20; -0.5, 5];
    z=[(1-e)*ones(1,20),ones(1,20),-0.5*[1,1,1,1,1]];
    p = poly(leja(z));
    
    z = [1-e, 20; 1, 20; -0.5, 5];
    
% JWF
%    if norm(imag(z(:,1))) == 0
%        fprintf('                 roots         multiplicities\n');
%        fprintf('\n');
%        fprintf('%25.15f \t \t \t %3g \n', z');
%    else
%        fprintf('                 roots ')
%        fprintf('   \t\t\t\t\t\t     multiplicities\n');
%        fprintf('\n');
%        fprintf('%22.15f + %22.15f i \t \t %3g \n', ...
%            [real(z(:,1)),imag(z(:,1)),z(:,2)]');
%    end;
    fprintf('It is suggested that multRoot use input gamma = 10.\n\n');
    z = displayRoots(z,2); % This also clips small parts and sorts

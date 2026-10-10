function [p,z] = bt02
%
% Brugnano and Trigiante
%
    p = poly([ones(1,10),-1,-1,i,-i,2]);
    z = [1,10; -1, 2; i, 1; -i, 1; 2, 1];
    
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
    z = displayRoots(z,2); % This also clips small parts and sorts
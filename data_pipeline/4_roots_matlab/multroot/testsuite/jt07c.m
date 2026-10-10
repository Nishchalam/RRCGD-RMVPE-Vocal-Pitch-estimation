function [p,z] = jt07c
%
%  test polynomial suggested by Jenkins and Traub
%
a = 10^(-8);
% JWF	
%  p = poly([.001, .01, .1, .1+a*i, .1-a*i, 1, -10]);
   z = [.001, .01, .1, .1+a*i, .1-a*i, 1, -10];
   p = poly(z);
   z = [z',ones(7,1)];
   y = [-10, 1; 1, 1; 0.01, 1; 0.001, 1; 0.1, 3];
   fprintf('\n');
   fprintf(' Ill-conditioned polynomial \n');
   fprintf(' It is constructed using\n\n')
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
   fprintf('\n');
   fprintf(' However, the polynomial is closer to having \n');
   fprintf('                 roots         multiplicities\n');
   fprintf('\n');
   fprintf('%25.15f \t \t \t %3g \n', y');
   
   
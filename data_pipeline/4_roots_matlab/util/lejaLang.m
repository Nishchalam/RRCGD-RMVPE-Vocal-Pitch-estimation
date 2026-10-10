% function [x_out] = leja(x_in)
% Marcus Lang's leja sort routine together with a bug fix.
%
% function [x_out] = leja(x_in)
%
%    Input:   x_in
%
%    Output:  x_out
%
%    Program orders the values x_in (supposed to be the roots of a 
%    polynomial) in this way that computing the polynomial coefficients
%    by using the m-file poly yields numerically accurate results.
%    Try, e.g., 
%               z=exp(j*(1:100)*2*pi/100);
%    and compute  
%               p1 = poly(z);
%               p2 = poly(leja(z));
%    which both should lead to the polynomial x^100-1. You will be
%    surprised!
%
%

%File Name: leja.m
%Last Modification Date: %G%	%U%
%Current Version: %M%	%I%
%File Creation Date: Mon Nov  8 09:53:56 1993
%Author: Markus Lang  <lang@dsp.rice.edu>
%
%Copyright: All software, documentation, and related files in this distribution
%           are Copyright (c) 1993  Rice University
%
%Permission is granted for use and non-profit distribution providing that this
%notice be clearly maintained. The right to distribute any portion for profit
%or as part of any commercial product is specifically reserved for the author.
%
%Change History:
% James W. Fox 10/8/03
% The original code reported that leja(1:5) was [5 2 4 1 3]. This is clearly wrong because the second
% value in any leja sorted array is the value farthest from the first value. In this case the
% second value must be 1. The reason this happened was because the orignal code computed the
% product of the distances as prod(a(1:l,l:n)). Unfortunately the first row contains the absolute
% values of the original values. I.e. the product he actually got was the product he really wanted
% times the absolute value of the point in question. There are two ways to fix this. Either set
% a(1,:)=1; before the "for" loop or change prod(a(1:l,l:n)) to prod(a(2:l,l:n),1). The latter option
% starts the multiplications at row 2. Note, you need to specify a second argument to prod.

% COMMENT: If you use this to leja sort rand(1,500), all products will
% underflow. You are multiplying hundreds of numbers together and all are
% <1.

function [x_out] = leja(x_in)
x = x_in(:).'; n = length(x);

a = x(ones(1,n+1),:);
a(1,:) = abs(a(1,:));

[dum1,ind] = max(a(1,1:n));  
if ind~=1
  dum2 = a(:,1);  a(:,1) = a(:,ind);  a(:,ind) = dum2;
end
x_out(1) = a(n,1);
a(2,2:n) = abs(a(2,2:n)-x_out(1));

for l=2:n-1
  % [dum1,ind] = max(prod(a(1:l,l:n)));  ind = ind+l-1; % Lang's original line
  [dum1,ind] = max(prod(a(2:l,l:n),  1));  ind = ind+l-1; % Added by JWF

  if l~=ind
    dum2 = a(:,l);  a(:,l) = a(:,ind);  a(:,ind) = dum2;
  end
  x_out(l) = a(n,l);
  a(l+1,(l+1):n) = abs(a(l+1,(l+1):n)-x_out(l));
end
x_out = a(n+1,:);





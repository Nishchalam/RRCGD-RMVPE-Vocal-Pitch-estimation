function ans = flip(x);

% REVERSE THE ORDER OF THE COEFFICIENTS IN X.

% REQUIRED INPUT:
%   x = Any array.

% OUTPUT:
%   ans = The reversed coefficients.

% USES:
%   Nothing.

% AUTHOR:
%   James W. Fox 4/22/2001

% COMMENTS
% This is useful when applied to a complicated expression.
% For example, lroots has
%     z = [z, flip(conj(z(j)))];

ans = x(length(x) : -1 : 1);
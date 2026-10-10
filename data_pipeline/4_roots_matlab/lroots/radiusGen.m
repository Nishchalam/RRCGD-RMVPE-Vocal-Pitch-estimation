% function radii = radiusGen(coeff, radialExtent, radialM);
%
% PURPOSE:
%   Generate the radii, from 1 down to 0, to be used for grid search.
%
% REQUIRED INPUT:
%   coeff - The polynomial coefficients.
%   radialExtent - A value > 1. Increase this to increase the region where we
%                  have a tighter grid.
%   radialM - Radial multiplier. A value > 1. Increase this to decrease the
%             radial grid size.
%
% OUTPUT:
%   radii = An optimized set of search radii.

% USES:
%   Nothing.
%
% AUTHOR:
%   James W. Fox 4/23/2001

function radii = radiusGen(coeff, radialExtent, radialM);

% This would be a nice place to use a winding number function to
% try to determine the kind of distribution we are dealing with. I.e.,
% is it random coefficient or zero-phase Butterworth, etc.

% This would also be a nice place to use a winding number function to
% determine the minimum and maximum radius and adjust the returned array
% accordingly. Note, there must be TWO search radii < the minimum root
% radius and two search radii > the maximum root radius.

% distribution is the kind of radius distribution.
% Random coefficient is the only thing I have implemented.

degree = length(coeff) - 1;
distribution = 1; % Random coefficient.
radii = [];
radius = 1;
while(radius > 0)
  radii = [radii, radius];
  if(distribution == 1)
    NNDist = NNDistR(radius, degree, radialExtent);
  else
    error('radiusGen method still under construction.');
  end
  radius = radius - NNDist/(3*radialM);
end
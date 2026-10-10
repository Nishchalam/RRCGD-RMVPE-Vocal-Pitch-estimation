% function dist = NNDistR(radius, degree, radialExtent);
%
% PURPOSE:
%   Given the polynomial degree and a root radius, what is the minimum
%   distance to the nearest neighboring root?
%
% REQUIRED INPUT:
%   radius - An array of root radii. We seek the minimum distance to the nearest neighbor.
%   degree - The polynomial degree
%   radialExtent - Increasing this lengthens the region where
%                  we have a smaller radial grid. Usually equals 1,
%                  but the program may increase it to 2, 4, or 8,
%                  if the grid search did not find enough roots.
%
% OUTPUT:
%   dist = The minimum distance to the nearest neighbor for a root
%         of radius "radius" in a random coefficient polynomial of
%         degree "degree."
%
% ALGORITHM
%   By looking at nearest neighbor plots of 25,000 upper half-plane roots
%   inside the unit circle for polynomials of degree 100, 1,000 and 10,000,
%   I determined that there are two regions of interest - near 1 and far from 1. 
%
%   For radii near 1.0 a lower bound was approximately
%     distance = .0122229*degree*(radius - 1 - 8.09017/degree)^2.
%
%   For radii far from 1.0 a good lower bound is approximately
%     distance = (1 - radius - 10/degree)
%   What is 'far from 1?' Empirical tests say the break point is
%   r = 1 - 18/degree. 18 is the break3 default in lroots(). However,
%   break3 is not passed to this routine. 18 is just hardwired here.

% USES:
%   Nothing.
%
% AUTHOR:
%   James W. Fox 4/22/2001
% 
% CHANGES
%     1/26/2006 Sid Burrus computed nearest neighbor data for degree 1e6 and 2e6.
%     This was used to generate completely new fit functions. The previous functions
%     only used nearest neighbor data for degree 100, 1,000, and 10,000.

function dist = NNDistR(radius, degree, radialExtent);

% I am not going to bother doing error checks on the arguments.

Lradius = length(radius);
dist = zeros(1,Lradius);
w1 = max(0,1 - radialExtent*3/degree); 
w2 = max(0,1 - radialExtent*18/degree);
w3 = max(0,1 - radialExtent*50/degree);
w4 = max(0,1 - radialExtent*100/degree);
w5 = max(0,1 - radialExtent*500/degree);
w6 = max(0,1 - radialExtent*1000/degree);
w7 = max(0,1 - radialExtent*5000/degree);
w8 = max(0,1 - radialExtent*10000/degree);

for k = 1:Lradius
  thisRadius = radius(k);
  if(thisRadius >= w1)
     window = 1;
  elseif(thisRadius >= w2)
     window = 2;
  elseif(thisRadius >= w3)
     window = 3;
  elseif(thisRadius >= w4)
     window = 4;
  elseif(thisRadius >= w5)
     window = 5;
  elseif(thisRadius >= w6)
     window = 6;
  elseif(thisRadius >= w7)
     window = 7;
  elseif(thisRadius >= w8)
     window = 8;
  else
     window = 9;
  end
  
  while(1) % The only reason for while is so I can have access to "break." This is only executed once.
    if(window == 1 | window == 2)
      x1 = w1; 
      if(degree <= 5000)
        y1 = 0.0710*log10(degree) - 0.4911;
      elseif(degree <= 10000)
        y1 = 0.0166*log10(degree) - 0.2914;
      else
        y1 = .0379*log10(degree) - .3834;
      end
      
      if(window == 1)
        y2 = y1;
        x2 = 1;
        break;   
      end
    end
   
    if(window == 2 | window == 3)
      x2 = w2;       
      if(degree <= 50)
        y2 = 2.3253*log10(degree) - 5.5507;
      elseif(degree <= 500)
        y2 = 0.5007*log10(degree) - 1.9015;
      elseif(degree <= 1000)
        y2 = 0.3322*log10(degree) - 1.4466;
      elseif(degree <= 10000)
        y2 = 0.1488*log10(degree) - 0.8973;
      else
        y2 = 0.0554*log10(degree) - 0.5183;
      end
      
      if(window == 2)
        break;
      end
    end
    
    if(window == 3 | window == 4)
      x1 = w3;
      if(degree <= 100)
        y1 = 11.62*log10(degree) - 25.75;
      elseif(degree <= 500)
        y1 = 2.4035*log10(degree) - 7.3071;
      elseif(degree <= 10000)
        y1 = 0.3345*log10(degree) - 1.6893;
      elseif(degree <= 100000)
        y1 = .1024*log10(degree) - .8037;
      else
        y1 = 0.0707*log10(degree) - 0.6332;
      end
         
      if(window == 3)
         break;
      end
    end
    
    if(window == 4 | window == 5)
      x2 = w4;
      if(degree <= 500)
        y2 = 5.3650*log10(degree) - 15.7301;
      elseif(degree <= 1000)
        y2 = 1.2291*log10(degree) - 4.5673;
      elseif(degree <= 10000)
        y2 =  0.4437*log10(degree) - 2.2008;
      else
        y2 = 0.0764*log10(degree) - 0.6884;
      end
        
      if(window == 4)
         break;
      end
    end
         
    if(window == 5 | window == 6)
      x1 = w5;
      if(degree <= 1000)
         y1 = 7.6404*log10(degree) -25.4213;  
      elseif(degree <= 5000)    
         y1 = 2.3034*log10(degree)  -9.4102;
      elseif(degree <= 10000)
         y1 = 0.6312*log10(degree) -3.2247;
      elseif(degree <= 100000)
         y1 = 0.3000*log10(degree) - 1.900;
      else
         y1 = 0.1139.*log10(degree) - .9708;
      end
      
      if(window == 5)
         break;
      end
    end
         
    if(window == 6 | window == 7)
      x2 = w6;
      if(degree <= 5000)    
         y2 = 4.9788*log10(degree) - 19.7363;
      elseif(degree <= 10000)
         y2 = 1.3454*log10(degree) - 6.2965;
      elseif(degree <= 100000)
         y2 = 0.4350*log10(degree) - 2.6550;
      else
         y2 = 0.1497*log10(degree) - 1.2087;
      end
      
      if(window == 6)
         break;
      end
    end
         
    if(window == 7 | window == 8)
      x1 = w7;
      if(degree <= 10000)
         y1 = 8.6370*log10(degree) - 36.9481;  
      elseif(degree <= 100000)
         y1 = 1.7000*log10(degree) - 9.2000;   
      else
         y1 = 0.2650*log10(degree) - 2.0044;
      end
      
      if(window == 7)
         break;
      end
    end
    
    if(window == 8 | window == 9)
      x2 = w8;
      if(degree <= 500000)
         y2 = 0.4494*log10(degree) - 3.1259 ;   
      elseif(degree <= 1e6)
         y2 = 0.3089*log10(degree) - 2.3256 ;
      else
         y2 = 0.1727*log10(degree) -1.5084 ;
      end
      
      if(window == 8)
         break;
      end
    end
    
    if(window == 9)
      x1 = 0;
      y1 = -5.0;
      break;
    end    
  end
  
  % Find y=mx+b that goes through (x2,y2) & (x1,y1)
   m = (y1 - y2) / (x1 - x2);
   b = y1 - m*x1;
  dist(k) = 10^(1./(m*thisRadius + b));
end

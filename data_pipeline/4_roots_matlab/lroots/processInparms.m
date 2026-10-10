% function [ID_NUM, radialExtent, thetaMult, radialMult, break1, break2, break3] = processInparms(inParms)
%
% PURPOSE:
%   Set defaults for inParms and see if the user wants to over-ride some of them.   
%
% OPTIONAL INPUT:
%   inParms - The (optional) second argument to lroots.
%
% OUTPUT:
%   See lroots argument documentation.   

% USES:
%   Nothing.
%
% AUTHOR:
%   James W. Fox 12/28/05

function [ID_NUM, radialExtent, thetaMult, radialMult, break1, break2, break3] = processInparms(inParms)

%---- Set defaults for inParms ----
ID_NUM = 1 ;
radialExtent = 1;
radialMult = 1 ;
thetaMult = 1 ;
break1 = 5;   % Breaks where thetaM is decreased by a factor of 2.
break2 = 10;  % This is in units of inc=1/degree. I.e. break1=5*inc.
break3 = 18;

if(nargin) % User wants to over-ride some of the defaults.
  ID_NUM = inParms(1);
  NinParms = length(inParms);

  if(NinParms > 1)
    radialExtent = inParms(2);
    if(radialExtent <= 0)
      error('$2(2) must be > 01.');
    end
  end

  if(NinParms > 2)
    thetaMult = inParms(3);
    if(thetaMult <= 0)
      error('$2(3) must be > 0.');
    elseif(thetaMult ~= 2^nextpow2(thetaMult))
      error('$2(3) must be a power of 2.');
    end
  end

  if(NinParms > 3)
    radialMult = inParms(4);
    if(radialMult <= 0)
      error('$2(4) must be > 0.');
    end
  end

  if(NinParms > 4)
    break1 = inParms(5);
    if(break1 ~= round(break1) | break1 < 1)
      error('$2(5) must be a postive integer.')
    end
  end

  if(NinParms > 5)
    break2 = inParms(6);
    if(break2 ~= round(break2) | break2 < 1)
      error('$2(6) must be a postive integer.')
    end
  end

  if(NinParms > 6)
    break3 = inParms(7);
    if(break3 ~= round(break3) | break3 < 1)
      error('$2(7) must be a postive integer.')
    end
  end

  if(break2 < break1)
    error('Must have $2[5] <  $2[6]');
  elseif(break3 < break2)
    error('Must have $2[6] < $2[7]');
  end
end

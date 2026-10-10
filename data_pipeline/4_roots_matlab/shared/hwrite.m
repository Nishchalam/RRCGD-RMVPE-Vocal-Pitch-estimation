% function hwrite(s);
%
% PURPOSE:
%   Write a character string to the history log, roots.his
%
% REQUIRED INPUT:
%   s = Any character string.
%
% OUTPUT:
%   None.

% USES:
%   global HISTORYFILE.
%
% AUTHOR:
%   James W. Fox 4/22/2001

function hwrite(s);

global HISTORYFILE; % Set in rootsIni

fp = fopen(HISTORYFILE,'a');
if(fp)
  fprintf(fp,'\n%s',s);
  fclose(fp);
end
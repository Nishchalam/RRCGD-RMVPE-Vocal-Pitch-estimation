% function [error] = tunwrap2(deg, seeds )
%
% PURPOSE
%    Return the error bewteen a single call to unwrapm and phase unwrapping by factorization.

% AUTHOR
% James W. Fox 2/8/2003

function error = tunwrap2(deg, seeds)
%fp=fopen('tunwrap.txt','w');
%if(fp)
   Nseeds = length(seeds);
   error=zeros(1,Nseeds);
   for k=1:Nseeds
      fprintf(1,'\n\t\t\t<<<<<<Starting iteration %d>>>>>>\n',seeds(k));
      rand('seed', seeds(k));
      x=rand(1,deg+1);
	   [z,zerr]=lroots(x);
      [ph, amp, err]=unwrapz(x,z,zerr);
      p = unwrapm(x);
      error(k) = max(abs(p-ph));
%      fprintf(fp,'Seed=%d \t Pad=%d \t Error=%g pi\n',seeds(k), error(k));
%   end
%   fclose(fp);
%else
%   error('Could not open output file tunwrap.txt');
end

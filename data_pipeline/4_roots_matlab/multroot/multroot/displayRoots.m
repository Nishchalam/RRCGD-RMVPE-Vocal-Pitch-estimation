function z2d = displayRoots(z2d,flag)

% REQUIRED INPUT ARGUMENTS
%  z2d - a 2D array - roots in column 1 and multiplicities in column 2.
%
% OPTIONAL INPUT ARGUMENTS
%  flag - Default=1. 
%         0 - Do not do any printing, just clip small parts and sort.
%         1 - Sort and print with column heading "computed roots."
%         2 - Sort and print with column heading "actual roots."
%             "Actual roots" should be used if this is called from testsuite programs.

   if(nargin == 0)
      error('Supply 1 or 2 arguments.');
   elseif(nargin == 1)
      flag=1;
   elseif(flag ~=0 & flag ~=1 & flag~=2)
      error('$2 must equal 0, 1 or 2');
   end
   
   % For example, in uhlig04, clipSmallParts turns
   % 2.727272727272728e-001 + 2.160673918824248e-016 i    into
   % 2.727272727272727e-001 + 0.000000000000000e+000 i.
   z2d = clipSmallParts(z2d); % This is critically important before sorting.
   % First sort by increasing argument. Then sort positive/negative purely real/imaginary
   % values by magnitude. Do the same sort on the second column.
   [junk,j] = csort(z2d(:,1),2);
   z2d(:,1) = z2d(j,1);
   z2d(:,2) = z2d(j,2);
   if(flag == 0)
      return;
   end

   % BATCH MODE: disabled diagnostic output
   % rts = z2d(:,1);
   % multiplicities = z2d(:,2);
   % FID = 1; % stdout; i.e. the screen
   % if(~any(imag(rts))) % All roots are real.
   %    if(flag == 1)
   %       fprintf(FID,'        computed roots         multiplicities\n');
   %    else
   %       fprintf(FID,'        actual roots         multiplicities\n');
   %    end
   %    fprintf(FID,'\n');
   %    fprintf(FID,'%25.15e             %3g \n', z2d');
   % else
   %    if(flag == 1)
   %       fprintf(FID,'        computed roots ');
   %    elseif(flag==2)
   %       fprintf(FID,'        actual roots ');
   %    end
   %    fprintf(FID,'                              multiplicities\n');
   %    fprintf(FID,'\n');
   %    fprintf(FID,'%22.15e + %22.15e i \t %3g \n', ...
   %           [real(rts),imag(rts), multiplicities]');
   % end
   % if(FID > 2)
   %    fclose(FID);
   % end
% function regtest
%
% PURPOSE
%   Multroot regression testing. The output has two printed columns. Column 1 is the name
%   of the programs being tested. Each program tested returns both the polynomial and 
%   the exact roots. Column 2 is the maximum difference between the exact roots and the
%   roots produced by multroot. If multroot produces a different number of unique roots 
%   (eg. it fails to find multiple roots when there are multiple roots) then the 
%   difference is Inf.
%   As of June 2004, there are 110 programs to test.
%
% INPUT ARGUMENTS
%   None
%
% AUTHOR
%   James W. Fox May-June 2004.

% The following line in the program determines what programs to test.
%     progNum = 1:nNames; % Indicies of programs to do regression testing on.
% This will test every program. If you only want to test fib05, which is the 
% sixth program, then change this to
% progNum = 6; 
% If you want to test programs 4,7, and 12 then change this to
% progNum = [4,7,12];
% If you want to test the last program then change this to
% progNum = nNames;

function regtest

% names is a list of all the files in the testsuite directory. Four are helper functions
% which do not return a coefficient array and a root array - fib, fl, lgd, and triple.
% These names are commented out.
names = [
'bt01;     ';
'bt02;     ';
'bt03;     ';
'bt04;     ';
'farloi01; ';
'fib05;    ';
'fib10;    ';
'fib15;    ';
'fib20;    ';
'fib30;    ';
'fib50;    ';
'fib100;   ';
'fib150;   ';
%'fib;      '; % Does not return a root array.
'fibsq04;  ';
'fibsq08;  ';
'fibsq16;  ';
'fibsq24;  ';
'fibsq32;  ';
'fibsq48;  ';
'fl01;     ';
'fl02;     ';
'fl03;     ';
'fl04;     ';
'fl05;     ';
'fl06;     '; % Fails about 1 time out of 7.
'fl07;     '; % Fails about 1 time out of 2.
%'fl;       '; % Does not return a root array.
'henrici;  ';
'igyp00;   ';
'igyp01;   ';
'igyp02a;  ';
'igyp02b;  ';
'igyp02c;  ';
'igyp03;   ';
'iliev00;  ';
'iliev01;  ';
'iliev02;  ';
'iliev03;  ';
'inex01;   '; % inex** round the coefficients to 4-7 significant digits
'inex02;   '; % and this destroys the multiple roots.
'inex03;   ';
'inex04;   ';
'jt01a;    ';
'jt01b;    ';
'jt02;     '; % multroot is much worse than roots; even finding a (false) double root.
'jt03;     ';
'jt04;     ';
'jt05;     ';
'jt06;     '; % multroot finds the wrong number of multiple roots.
'jt07a;    '; % multroot finds the wrong number of multiple roots.
'jt07b;    '; % multroot finds the wrong number of multiple roots.
'jt07c;    '; % multroot finds the wrong number of multiple roots.
'jt07d;    '; % multroot finds the wrong number of multiple roots.
'jt08;     ';
'jt09;     ';
'jt10a;    ';
'jt10b;    ';
'jt10c;    ';
'jt11a;    '; % Roots MUST be sorted by leja before unfactoring or results are poor.
'jt11b;    '; % Roots MUST be sorted by leja before unfactoring or results are poor.
'jt11c;    '; % Roots MUST be sorted by leja before unfactoring or results are poor.
'large01;  '; % Coefficients are intentionally rounded up at the 10-th digit. Results are ok.
'large02;  '; % Coefficients are intentionally rounded up at the 10-th digit. Results are ok.
'large03;  '; % Coefficients are intentionally rounded up at the 10-th digit. Results are ok.
'large04;  '; % Coefficients are intentionally rounded up at the 10-th digit. Results are poor.
'large05;  '; % Coefficients are intentionally rounded up at the 10-th digit. Results are poor.
'lgd05;    ';
'lgd10;    ';
'lgd15;    ';
'lgd20;    ';
'lgd24;    ';
'lgd50;    '; % Forward error = 167! Max error = 2.00
'lgd100;   '; % Forward error = 77200! Max error = 2.22
%'lgd;      '; % Does not return a root array.
'miyak00;  ';
'miyak02;  ';
'miyak04;  ';
'miyak08;  ';
'near01;   ';
'near02;   ';
'near03;   '; % This fails about 20% of the time.
'petk01;   ';
'petk02;   ';
'petk03;   ';
'petk04;   ';
'petk05;   ';
'petk06;   ';
'petk07;   ';
'spiral;   ';
'toh01;    ';
'toh02;    ';
'toh03;    '; % Extremely ill-conditioned. Roots are supposed to be real, but multroot says no.
'toh04;    ';
'toh05;    ';
'toh06a;   ';
'toh06b;   ';
'toh06c;   ';
'triple01; ';
'triple02; ';
'triple03; '; % Fails about 1 time out of 2.
'triple04; '; % Fails about 1 time out of 2.
%'triple;   '; % Does not return a root array.
'twin01;   ';
'twin02;   ';
'twin03;   ';
'twin04;   ';
'uhlig01;  ';
'uhlig02;  ';
'uhlig03;  ';
'uhlig04;  ';
'uhlig05;  ';
'uhlig06;  ';
'underflow;';
];

[nNames,nameLength]=size(names);
% progNum is the indicies of programs to do regression testing on.
progNum = 1:nNames; % Indicies of all programs.
% Uncomment the next line if you only want to test the difficult polynomials.
% progNum = [25 26 44 48:52 64:65 70:71 88 97]; % Indicies of difficult polynomials.
% inex01-inex04 are inexact polynomials; they round their coefficients to 4-10
% significant digits. Of course this destroys multiple roots. There is no way to get
% accurate roots if you do this. I am going to ignore them and not count them as 
% "difficult polynomials." The above list is the other "difficult polynomials" that should
% work but do not. It would be worthwhile to play with multroot's the last 3 optional
% input arguments to see if changing them can make these polynomials work.
err = zeros(1,progNum(end)-progNum(1)+1); % Allocate space for the error array.

for k=progNum
   n = deblank(names(k,:));
   fprintf(['\n--------------\n', n,'\n']);
   [p,z]=feval(n);
   [pp,zz]=multroot(p);
   % The new version of multroot adds an output argument. Test if we are using the
   % old multroot and make adjustments to the output. This allows it to run but the
   % old multroot code does not call displayRoots so the roots are in the wrong order.
   if(length(p)>2 & length(zz)==1) % We are using the old code.
      zz = pp;
   end
   if(length(z) ~= length(zz)) % multroot found the wrong number of distinct roots.
      err(k) = Inf;
   else
      err(k) = max(max(abs(z-zz)));
   end
end

fprintf('\nFunction name   Maximum multroot error\n');
fprintf('-------------   ---------------------\n');
for k=progNum
   n = deblank(names(k,:));
   if(err(k) < 1e-8)
      fprintf('%s \t%g\n',names(k,1:end),err(k));
   else
      fprintf('%s \t========> %g\n',names(k,1:end),err(k));
   end
end

m = find(err >= 1e-8);
if(~isempty(m))
   fprintf('\nThe worst programs and their errors were:\n');
   fprintf('\nFunction name   Maximum multroot error\n');
   fprintf('-------------   ---------------------\n');
   for n=1:length(m)
      fprintf('%s \t%g\n',names(m(n),1:end),err(m(n)));
   end
end

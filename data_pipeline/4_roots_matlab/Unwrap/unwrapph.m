% function ph = unwrapph(x,N);
%
% PURPOSE
%   Phase unwrapping by zero padding once with a lot of zeros. Compare this
%   with unwrapm which starts out with a small zero padding, then repeatedly
%   doubles the padding, and stops when two consecutive answers agree at the
%   central frequency.
%
%   For 1,000-degree random coefficient polynomials, if you use the defaults,
%   unwrapm gets the correct answer, unwrapz's answer, 85% of the time. 
%   Unwrapph gets the correct answer 98.5% of the time.
%
%   If there is a root on the unit circle, only unwrapz can get the correct answer.
%
% REQUIRED INPUT ARGUMENTS
%  x - Any time series
%
% OPTIONAL INPUT ARGUMENTS
%  N - Default = 128 if degree<=1000
%              = 256 otherwise. 
%      Pad with zeros until the length is N*length(x). This seems to get the same
%      answer as unwrapz at least 98.5% of the time if degree<4096. Higher degrees
%      require larger value of N. Lower degrees could use smaller values.
%
% REQUIRED OUTPUT ARGUMENTS
%  ph - The unwrapped phase
% 
% ALGORITHM
%  Traditional unwrappers can fail because the actual phase jumps can be larger than pi.
%  This can be remedied by first padding with zeros before taking the FFT. This
%  subsamples the phase at a smaller step size so the phase jumps are smaller and
%  satisfy the assumptions of traditional unwrappers.
%
% USES
%    angl, adjustUnwrapping

% AUTHOR:
%    James W. Fox

function ph = unwrapph(x,N);

if(nargin < 1)
   error('Supply $1.');
elseif(min(size(x)) ~= 1)
   error('$1 must be 1xn or nx1.');
elseif(nargin < 2 | isempty(N))
   if(length(x) <= 1001)
      N = 128;
   else 
      N = 256;
   end
elseif(~isreal(N) | N<1 | round(N) ~= N | max(size(N)) ~= 1)
   error('$2 must be a real positive integer.');
end

xx = [x(:)', zeros(1, (N-1)*length(x))];
ph = unwrap(angl(fft(xx)));
ph = reshape(ph(1:N:length(ph)),  size(x));
s = fft(x);
ph = adjustUnwrapping(ph, abs(s), angl(s));

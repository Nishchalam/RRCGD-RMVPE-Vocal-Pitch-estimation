% function [ph, amp, err] = unwrapz(x, z, zerr, N)
%
% PURPOSE
%   Use a real polynomial's roots to unwrap its amplitude and phase spectra.
%   To inverse transform: y = ifft(amp.*complex(cos(ph),sin(ph)));
%
%   If 1 is not a root, this function returns the one unwrapped phase
%   function that begins at 0. If there are no roots on the unit circle, 
%   the value at the central frequency is -pi*#(roots outside the unit circle)
%   which is connected to Gauss's Theorem on the winding number. 
%
%   For a given polynomial, there are infinitely many candidates for the unwrapped 
%   phase. If ph is any unwrapped phase function, you can add any integer 
%   multiple of 2*pi. Also, you can add any odd integer multiple of pi to the 
%   phase and multiply the amplitude spectrum by -1. 
%
% REQUIRED INPUT
%   x - The polynomial coefficients. First term is the highest order term.
%   z - The roots of x.
%
% OPTIONAL INPUT
%   zerr - Error estimates for each root. If either [] or no value is supplied
%          then zerr = 1e-11 + zeros(size(z));
%   N - Pad x with zeros on the right to length N. Why? Suppose y3 = conv(y1,y2)
%       and you wish to verify that the unwapped phase(y3) equals the sum of the
%       unwrapped phases of y1 and y2. You need to ensure the 3 unwrappings
%       have the same length.
%
% OUTPUT
%   ph -  The unwrapped phase spectrum.
%         ph(1) == -pi/2*(the multiplicity, modulo 2, of 1.0 as a root).
%
%         ph(end)/(2*pi) is approximately
%         -#(roots outside the unit circle)
%         - 0.5*#(roots on the unit circle including 1.0) 
%         - 0.5*(the multiplicity, modulo 2, of 1.0 as a root)
%
%         It is not exactly equal because the FFT evaluates at the n'th roots
%         of unity and the last root of unity is close to but not equal to 1.
%         However, a related fact is true at the central frequency.
%         If n=length(ph) is even, then ph(1+n/2)/pi should equal this value. 
%         If n is odd, then mean(ph((n+1)/2:1+(n+1)/2)/pi should equal this
%         value. This fact is the basis for function unwrapm, an alternate
%         phase unwrapper that does not need to factor the polynomial.
%
%   amp - The unwrapped amplitude spectrum. Let z be any complex number with
%         z = r*exp(i*theta) where r and theta are real and r>0. Note that
%         z = -r*exp(i*(theta+pi)). There is no obvious reason for disallowing
%         this as a polar coordinate representation just because -r<0. In other
%         words, why not allow the amplitude spectrum to be sometimes/always
%         negative? Why is the useful? The primary benefit is that if there
%         are roots exactly on the unit circle, then the phase undergoes a
%         discontinuity of exactly pi. If you subtract pi (not 2*pi) from some
%         of the angles and multiply the usual amplitude spectrum by -1 for the
%         same frequencies, then you can remove this discontinuity.
%
%         Consider the following example with roots i and -i on the unit circle.
%         If you subtract pi from certain angles and 2*pi from other angles, the
%         result is essentially a perfectly straight line. This is surely a 
%         natural representation. Furthermore, the unwrapped amplitude spectrum
%         looks much more "natural" than the wrapped amplitude spectrum.
%
%         y=[1 0 1 zeros(1,50)];
%         s=fft(y);
%         amp=abs(s);
%         plot(amp) % Note the 2 non-differentiable points.
%         title('Wrapped Amplitude'); 
%
%         ph=angle(s);
%         plot(ph) % Note jumps of size pi.
%         title('Wrapped Phase');
%
%         amp(15:40) = -1*amp(15:40); % Unwrap the amplitude. <==
%         plot(amp) % Now it looks differentiable.
%         title('Unrapped Amplitude');
%
%         ph(15:40)=ph(15:40)-pi;   % Unwrap part of the phase by subtracting pi. <==
%         ph(41:53)=ph(41:53)-2*pi; % Unwrap the rest of the phase by subtracting 2*pi.
%         plot(ph);
%         title('Unrapped Phase');
%
%         plot(ph); This is indeed a straight line; the slope is constant.
%         format long e
%         [min(diff(ph)), max(diff(ph))]
%         	ans = -1.185506661732685e-001   -1.185506661731601e-001
%
%         Check that the unwrapped spectra reconstructs the original polynomial.
%         x=ifft(amp.*complex(cos(ph),sin(ph)));
%         max(abs(x-y))
%        	ans =  2.0892e-014

%   err - An indication of the maximum phase alteration made at any frequency so that the 
%         unwrapped phase obtained from the roots agrees with the traditional wrapped 
%         phase mod pi. 0<=err<=1.
%
%         0 - Complete success. Not even one bit of round-off error in any phase value.
%             For 100-degree seismic data, expect values like 1E-13. For 1,000-degree
%             seismic data, expect values like 1E-10.
%         1 - Total failure. Anything close to 1 indicates that the original unwrapping
%             was extremely poor.
%
%             A 100 degree polynomial was factored to get a complete set of roots, z.
%             Unwrapping with z gave err==9E-14. This was the correct unwrapped phase.
%
%             Unwrapping with z+.0001 gave err==.015 but the phase was identical with
%             the phase obtained using z. In other words,the correct unwrapped phase 
%             was obtained even though every root had an error of .0001.
%
%             Unwrapping with z+.001 gave err=.18. The unwrapped phase was incorrect.
%
%         First a tentative unwrapped phase is obtained from the roots. The wrapped
%         phase is assumed to be perfectly correct modulo pi. Thus at each frequency
%         some value in the range -pi/2 to pi/2 is added to the tentative unwrapped
%         phase to force it to agree with the wrapped phase mod pi. Thus, if the
%         unwrapped phase obtained from roots has absolute errors less than pi/2
%         this procedure automatically corrects those errors. Thus the correct answer
%         can often be obtained even if some roots are ill-conditioned and can not be 
%         obtained with great precision. That is how the correct answer was obtained
%         above when z+.0001 was unwrapped instead of z.
%
%         Why should the wrapped and unwrapped phase agree mod pi and not
%         mod 2*pi? The polynomials x and -x have the same set of roots. Adding pi to
%         every value of the phase multiplies the time series by -1. Thus you could
%         argue that unwrapped phases of x and -x differ by pi. However, since I allow
%         negative values in the amplitude spectrum, let us take the approach
%         that x and -x have identical unwrapped phase functions but their unwrapped
%         amplitude specra differ by a factor of -1. One nice thing about this is that
%			 all unwrapped phase functions begin at 0 (unless 1 is a root).
%
% SAMPLE USAGE
% 		x = rand(1,101);
%		[z,e] = lroots(x);
%		[ph,amp,err] = unwrapz(x,z,e);
% 		plot(ph)

% USES
%    angl, ixFirstNZ, adjustUnwrapping

% AUTHOR
%   James W. Fox. 2/15/02

% COMMENTS
%   ZERO PADDING:
%     Zero padding on the right is innocuous for siesmic data but not for Matlab
%     polynomials.

%		For Matlab polynomials, padding zeros on the left of an array does not
%     change the polynomial because, in Matlab conventions, the first coefficient
%     is the highest order term. However padding zeros on the right changes the
%     polynomial degree.

%     On the other hand, padding zeros on the left of a seismic trace changes
%     the time of every event on the trace. Only padding on the right is natural.

%     In Numerical Recipes' conventions, the first coefficient is the constant term.
%     Thus, like seismic data, zero padding on the right is innocuous. For this reason
%     and others listed below, Numerical Recipes' conventions are easier than Matlab's
%     conventions for the purposes of phase unwrapping.

%     If you have a 1,000 point seismic trace and you pad 24 points on the right, you
%     would expect the two traces to have the same unwrapped phase functions. This
%     will indeed be true of the algorithm below. However, this also means that
%     the algorithm below will give the same unwrapped phase function to
%     [1 2 3 4 5 6] and [1 2 3 4 5 6 0 0] although they are quite different when
%     viewed as Matlab polynomials. The latter has 0 as a root of multiplicity 2.

%   GAUSS'S THEOREM:
%     Gauss's Theorem says that a certain integral equals 2*pi*i*#(maxphase roots)
%     You might think the last value of the unwrapped phase would be the same but,
%     if there are no unit circle roots, the last value of the unwrapped phase is
%     approximately -2*pi*#(minphase roots).

%     Why "minphase" when Gauss has maxphase? Again, it is because of 
%     Matlab's conventions. If z is a maxphase root of x which is viewed as a Matlab
%     polynomial with first coefficient being the highest order coefficient, then
%     1/z is a minphase root of x viewed as a Numerical Recipes polynomial with first
%     term being the constant term. If we use Numerical Recipes' conventions
%     then it would have been maxphase, just like Gauss.

%     Why is it -2*pi while Gauss has +2*pi? It depends on how you do your FFT.
%     Some people use exp(-2*pi*i*omega*t) for the Fourier kernel and others use
%     exp(+2*pi*i*omega*t). Mathematicians and Matlab use the former while 
%     electrical engineers and Numerical Recipes uses the latter. Switching from
%     one version of the Fourier kernel to the other causes a multiplication by
%     -1 in phase unwrapping.

%		Where is the "i?" The Gauss integral is not being computed. Rather something
%     related to this integral is being computed.

%   LINDSEY'S TRICK
%     Joe Pat Lindsey discovered Steiglitz and Dickinson's article "Phase Unwrapping
%     by Factorization" and recognized that his work with factorization of polynomials
%     might also lead to a better way to unwrap phase. He further suggested that the phase
%     associated with a root r was computable from the negative of the phase associated
%     with 1/r. This would be useful because if the phase of r wrapped, the phase of 1/r
%     would not wrap and could easily be computed. Then the phase of 1/r could be used 
%     to compute the phase of r which wraps. In other we can easily compute the phase of
%     a root, r, that wraps from 1/r which does not wrap.
%
%     Fox worked out the details and also discovered how to handle roots on the unit
%     circle and amplitude spectrum unwrapping.

%		Consider a real root r and the associated polynomial factor z-r which is denoted
%     as f = [1, -r, 0, 0, ..., 0, 0]. You must pad all the linear and quadratic
%     factors with zeros out to the length of the original polynomial. Otherwise the
%     phase of [1,-r] would be a 2 term array.

%     By Gauss's Theorem, if |r|>1 the phase function will not wrap at all. In other 
%     words, angle(fft(f)) is the unwrapped phase. You do not need to do any additional
%     work to unwrap it.

%     Again, by Gauss's Theorem, if |r|<1 it wraps exactly once. Seemingly, you must do
%     some work to unwrap it. Lindsey's ideas show how to do it without any work.

%     Consider the polynomial (z-r)*(z-1/r) = z^2 -(r+1/r)*z + 1. Write this as
%     [1, -(r+1/r), 1, 0, 0, 0, ..., 0, 0]. If you circularly rotate this left one sample
%     you get [-(r+1/r), 1, 0, 0, ..., 0, 0, 1]. This function posseses even symmetry about
%     its first point. Thus its Fourier phase is identically zero. In other words, if we add
%     the linear phase function to z^2 -(r+1/r)*z + 1 that corresponds to left shifting one
%     sample, then we get zero phase. This linear phase function is 2*pi*(0:N-1)/N where
%     N=length(x). To see that this does indeed left shift by one sample:
%     >> x=1:8;
%     >> N=8;
%     >> a=abs(fft(x));
%     >> p=angle(fft(x));
%     >> pshift=2*pi*(0:N-1)/N;
%     >> ifft(a.*complex(cos(p+pshift),sin(p+pshift))) % [2 3 4 5 6 7 8 1]
%
%     ans =
%       Columns 1 through 4 
%        2.0000 + 0.0000i   3.0000 + 0.0000i   4.0000 - 0.0000i   5.0000 + 0.0000i
%       Columns 5 through 8 
%        6.0000 - 0.0000i   7.0000 - 0.0000i   8.0000 - 0.0000i   1.0000 + 0.0000i


%     Thus,
%     phase(z^2 -(r+1/r)*z + 1) + 2*pi*(0:N-1)/N = 0 phase. Equivalently,
%     phase(z^2 -(r+1/r)*z + 1) = -2*pi*(0:N-1)/N. Thus,
%		phase(z-r) + phase(z-1/r) = phase(z^2-(r+1/r)*z+1) = -2*pi*(0:N-1)/N.
%     Thus phase(z-r) = -phase(z-1/r) - 2*pi*(0:N-1)/N.

%     Thus we can compute the phase corresponding to 1/r, which does not wrap, 
%     and use it to compute the phase corresponding to r which wraps once. 
%     The one wrap is accounted for by subtracting the linear phase.

%     Now consider a general complex root r which is not on the unit circle. If the 
%     polnomial has real coefficients, both r and conj(r) are roots so we want to 
%     add the phase contribution of (z-r)*(z-conj(r)). By Gauss, if |r|>1 then 
%     phase((z-r)*(z-conj(r))) does not wrap at all. Its unwrapped phase is 
%     angle(fft((z-r)*(z-conj(r)))). Of course, you must zero pad before taking
%     the fft.

%     What if |r|<1? Consider the polynomial (z-r)*(z-conj(r))*(z-1/r)*(z-conj(1/r)). 
%     This is a fourth degree polynomial. It has 5 coefficients and is even symmetric 
%     about the center term. Thus if you add a linear phase function corresponding to a
%     left shift by 2 samples, the resulting polynomial is zero phase. The linear phase 
%     function is 4*pi*(0:N-1)/N. Thus, like the above,
%     	phase((z-r)*(z-conj(r))) = -phase((z-1/r)*(z-conj(1/r))) - 4*pi*(0:N-1)/N.
%     Everything on the right side is easy to compute because the phase of the first
%     polynomial on the right does not wrap. The phase on the left wraps twice, 
%     corresponding to the two roots, because of the linear phase function on the right.

%     The general strategy is to sum the phase contribution of each root by:
%        If |r|>1, depending on whether r is real or not, use angle(fft) or
%     angle(fft((z-r)*(z-conj(r))) to compute the unwrapped phase of z-r and add
%     this to the sum. The phase does not wrap.
%        If |r|<1, compute the unwrapped phase of z-1/r or (z-1/r)*(z-conj(1/r)),
%     which does not wrap, and subtract it from the sum. Also subtract a linear
%     phase function which causes a wrap.

%     There is a "mistake" in all the above. Because of Matlab's convention
%     that the first term is the highest order term, roots outside the
%     unit circle wrap, not roots inside as Gauss would have you believe.

%     What about roots on the unit circle? In each case linear or quadratic
%     polynomials were constructed with roots on the unit circle and this was
%     padded with zeros. Then it was observed what was necessary to unwrap
%     the phase. It was noted that in each case the phase could be "unwrapped"
%     to get a perfectly straight line. However "unwrapped" might mean that it
%     was necessary to subtract (not 2*pi) from some of the angles. If
%     you do this, you must also multiply the amplitude spectrum by -1 at the 
%     corrresponding frequencies.

%     The following three cases assert that if a finite time series is circularly
%     symmetric about its first point then it is zero phase. That needs some
%     clarification. If it is circularly symmetric, then the spectrum is real; i.e.,
%     except for roundoff errors, it contains only positive and negative numbers.
%     Traditionally positive numbers have phase 0 while negative numbers have phase
%     pi so at first glance it does not seem to be zero phase. However, as mentioned
%     above, if you allow the amplitude spectrum to contain negative values, then
%     you can subtract pi from any phase value. This allow us to say that circularly
%     symmetric implies zero phase. However, you can only say this if you allow the
%     amplitude spectrum to be unwrapped.

%     CASE 1: r==-1. The associated polynomial is z+1 which is denoted 
%     [1 1 0 0 0 ... 0]. If you circularly shift left by 1/2 sample you get a
%     polynomial which is even symmetric about its first term. I.e., it is zero
%     phase. The linear phase shift corresponding to 1/2 sample left shift is
%     pi*(0:N-1)/N. Thus phase(z+1)+pi*(0:N-1)/N = zero phase. Equivalently, 
%     phase(z+1)=-pi*(0:N-1)/N.

%     CASE 2: r==1. The associated polynomial is z-1 which is denoted
%     [1 -1 0 0 0 ... 0]. If you circularly shift left by 1/2 sample you get a
%     polynomial which is odd symmetric about its first sample. If we do a 90
%     degree phase rotation, we get even symmetry which is zero phase. Thus 
%     phase(z-1) + pi*(0:N-1)/N + pi/2 = zero phase. Equivalently.
%     phase(z-1) = -pi*(0:N-1)/N - pi/2.

%     CASE 3: r is complex. Both r and conj(r) are roots so we need to add the 
%     phase contribution of (z-r)*(z-conj(r)) = z^2 - 2*real(r) + 1.
%     Write this as [1, -2*real(r), 1, 0 0 0 ..., 0]. Circularly shift this left
%     one sample to get a function with even symmetry about its first sample.
%     I.e., it is zero phase. Thus phase((z-r)*(z-conj(r))) = -2*pi*(0:N-1)/N.

%     However, there is a problem in CASE 3:
%     >> x=[1 0 1 zeros(1,13)]; % (z-i)*(z+i)=[1 0 1]
%     >> a=abs(fft(x))
%     a =
%       Columns 1 through 7 
%         2.0000    1.8478    1.4142    0.7654         0    0.7654    1.4142
%       Columns 8 through 14 
%         1.8478    2.0000    1.8478    1.4142    0.7654         0    0.7654
%       Columns 15 through 16 
%         1.4142    1.8478
%
%     >> p=angle(fft(x))
%
%     p =
%       Columns 1 through 7
%         0   -0.3927   -0.7854   -1.1781         0    1.1781    0.7854
%       Columns 8 through 14 
%         0.3927         0   -0.3927   -0.7854   -1.1781         0    1.1781%
%       Columns 15 through 16 
%         0.7854    0.3927
%
%     >> s=-2*pi*(0:15)/16;
%     >> p-s
%
%     ans =
%       Columns 1 through 7 
%         0   -0.0000         0   -0.0000    1.5708    3.1416    3.1416
%       Columns 8 through 14 
%         3.1416    3.1416    3.1416    3.1416    3.1416    4.7124    6.28325
%       Columns 15 through 16 
%         6.2832    6.2832
%
%     The first 4 values of p and s agree but values 6:12 differ by pi and values
%     14:16 differ by 2*pi. Values 5 and 13 are special because the amplitude
%     spectrum is zero here so the phase is actually undefined. If we multiply
%     the amplitude spectrum by -1 for samples 6:12, then we can add or subtract
%     pi from the corresponding phase values. This will make p and s agree at these
%     samples. This is unwrapping the amplitude as well as the phase spectrum.
%
%     This problem with CASE 3 applies to any complex root on the unit circle, not just
%     to i and -i. The amplitude spectrum must be unwrapped, i.e., multiply parts
%     of the amplitude spectrum by -1. If this is done, then the phase contribution
%     of these roots becomes a perfectly linear function with no discontinuities.

%     Allowing negative values in the amplitude may seem wrong. Let me try
%     to convince you with an example. Consider 17 points equally spaced on the
%     line imag(z)=real(z) in the complex plane with -1 <= real(z) <= 1. The first
%     point is 1+i in the first quadrant and the 17'th point is -1-i in the third
%     quadrant. The 9'th point is 0+0*i. If this set of points, is called pts, you
%     can do ifft(pts). I.e. pts is the discrete Fourier spectrum of some function.
%     (I do not care what the function is.) We already have the complex spectrum
%     of this function; it is the 17 points described above. What are the amplitude
%     and phase spectra? 

%     First let's use the traditional definitions. The first value of the spectrum
%     is 1+i so the amplitude is sqrt(2). The amplitudes decrease linearly to the
%     9'th point where the amplitude is 0. For the last half of the spectrum, the
%     amplitude increases linearly back up to sqrt(2). In other words, the
%     amplitude spectrum has a perfect "v" shape. The phase is 45 degrees for the
%     first half of the spectrum, zero at the midpoint, and it is -135 degrees for
%     the last half of the spectrum. The amplitude is non-differentiable at its
%     center point while the phase is discontinuous at the center point.

%     Now consider allowing negative values in the amplitude spectrum. The amplitude
%     spectrum starts at sqrt(2). It descends linearly to 0 at the midpoint. Then it
%     goes negative and continues to descend linearly to -sqrt(2). The graph of the
%     amplitude spectrum is a perfectly straight line. The non-differentiable point
%     has been eliminated. Meanwhile the phase spectrum is a constant 45 degrees
%     at all frequencies. This eliminates the point of discontinuity. This seems
%     more natural than the usual definitions.

%   FALIURE IN ADJUSTING THE UNWRAPPED PHASE:
%     In discussing the output argument, err, it was mentioned that a 100-degree 
%     polynomial was factored  to get a set of roots z. When the phase was unwrapped
%     using z+.0001 instead of z, it still got the correct answer.

%     Do not assume that we can always recover from errors smaller than .0001.
%     Suppose there is an ill-conditioned root very near the unit circle. You polish
%     it once and it lands slightly inside the unit circle. You polish it again and 
%     it lands slightly outside the unit circle. If you supply error estimates for
%     the roots, the program is supposed to detect this and decide that the root is
%     exactly on the unit circle. This may or may not be true and it makes a big 
%     difference in unwrapped phase whether the root is inside/outside/on the unit circle.

%   THEORY VERSUS EXPERIENCE:
%		The above discussion may seem like everything was done based on theory.
%     The reverse was true. Everything was done based on experience and then
%     a theory propounded to justify the experience. If you do anything different,
%     you must check	with experience. For example if you change conventions to 
%     let the first coefficient be the constant term or if you change the Fourier
%     kernel, then you need to get some experience.

% 		By experience, I mean: take a root r, and its conjugate if it is complex, and
%     form the linear or quadratic factor. Pad with zeros, compute the Fourier
%     spectrum, and compute the phase. Look to see what, if anything, needs to be done
%     to unwrap the phase. You need to check 10 cases. 1) Positive real maxphase 
%     2) Negative real maxphase 3) Complex maxphase 4) Positive real minphase 
%     5) Negative real minphase 6) Complex minphase 7) -1 as a root 8) 1 as a root
%     9) A complex root on the unit circle 10) 0 as a root.

%		Theory said the phase of a complex root on the unit circle would have a linear
%     phase function. It was only experience that showed it was also necessary
%     to unwrap the amplitude spectrum. 

%		Theory, i.e., Gauss's theorem, suggested that the phase associated with
%     a minphase root would wrap once while the phase for a maxphase root would
%     not wrap. Experience showed that if Matlab conventions were used, it would
%     be the other way around.

%   CONVENTIONAL PHASE UNWRAPPERS:
%		Conventional phase unwrappers are several orders of magnitude faster than
%	   factoring a polynomial and using the algorithm below. Unfortunately, quickly
%     getting an incorrect answer is not the best outcome. However, perhaps
%     traditional unwrappers can be greatly improved. The above discussion noted that if
%     1.0 is not a root, the first unwrapped phase value should be 0 and the value at the
%     central frequency is computable from #(minphase roots), the number of roots on the
%     unit circle, and the multiplicity of 1. As long as there are no roots on the unit
%     circle, there are winding number functions that can compute #(maxphase roots)
%     without factoring. Thus the central frequency value of the unwrapped phase is
%     a constraint that traditional phase unwrappers do not take into account. Perhaps
%     if they did, their results would be closer to the correct answer; perhaps
%     entirely correct.

function [ph, amp, err] = unwrapz(x, z, zerr, N)

% A later section checks that length(z) == length(x)-1. I am using
% Matlab's convention that the first coefficient is the highest order
% term, so leading zeros do not affect the polynomial. However they
% will affect length(x)-1. Therefore I must remove leading zeros.

k = ixFirstNZ(x); % Index of the first non-zero value.
if(k > 1)         % Trim off leading zeros.
   x = x(k : length(x));
end

if(nargin < 2)
   error('Must have at least 2 input arguments.');
elseif(any(imag(x)))
   error('$1 must be a real 1D array.');
elseif(~any(x))   
   error('Polynomial is identically zero.');
elseif(isempty(z) | length(z) ~= length(x)-1)
   error('Must have length($2) = (length($1) - 1) > 0');
elseif(nargin < 4)
   N = length(x);
elseif(isempty(N) | length(N)>1 | N < length(x))
   error('$4 must be a single positive integer that is >= length($2).');
else
   N = round(N); % In case they used an expression that had roundoff error.
end

if(nargin <= 2 | isempty(zerr))
   zerr = 1e-11 + zeros(size(z));
elseif(length(zerr) ~= length(z))
   error('Must have length($2)==length($3)');
end

k = find(imag(z) >= 0);
z = z(k);                     % Use only upper half-plane roots.
zerr = zerr(k);
ISODD = bitand(N,1);          % Is N odd or even?

NCxMin = 0;    % The number of complex minphase roots (not on unit circle).
NUCR = 0;      % The number of unit circle roots.
N_1 = 0;       % The number of times 1 is a root.
N_Minus_1 = 0; % The number of times -1 is a root.

sizeX = size(x); % Used at the very end.
x = x(:).';    % This is needed by the next line.
x = [x, zeros(1,N-length(x))]; % Pad to length N.
spectrum   = fft(x);
amp        = abs(spectrum);  % The usual amplitude spectrum. It will be unwrapped if necessary.
wrapped_ph = angl(spectrum); % The usual wrapped phase function.

if(ISODD)
   Nspectrum = (N+1)/2 + 1; 
else
   Nspectrum = N/2 + 1;
end
spectrum = spectrum(1:Nspectrum); % Truncate to half the spectrum.
% The spectrum is symmetric about sample Nspectrum. It is a waste of time to
% unwrap the second half when we can unwrap the first half and use symmetry
% to generate the second half. Thus the linear phase functions constructed 
% below are not identical with the long discussion above. They are half as long.

ph = zeros(1,Nspectrum);

% Table of exp(-2*pi*i*j/Nspectrum) j = 0,1,2,...,Nspectrum-1
theta = [0:(Nspectrum-1)]*(2*pi/N);
trigTable = complex(cos(theta), -sin(theta)); % exp(-i*theta)

% Matlab has tiny errors in their trig functions. Eg. sin(pi)=-1.2246e-16*i
% and cos(pi/2)=6.1232e-017. Correct these errors in trigTable
if(~mod(N,4)) % if(N is a multiple of 4)
   trigTable(N/4 + 1) = -i;
end
   
if(~mod(N,2)) % if(N is a multiple of 2)
   trigTable(Nspectrum) = -1;
end

for k = 1:length(z)
   z1 = z(k);
   absZ1 = abs(z1);
   imagZ1 = imag(z1);
   err1 = 2*zerr(k);
   
   if(abs(1 - absZ1) <= err1) % z1 is on the unit circle.
      if(~imagZ1) % z1 == 1 or -1
         if(abs(1 + z1) < err1) % z1 == -1
            NUCR = NUCR + 1;
            N_Minus_1 = N_Minus_1 + 1;
         elseif(abs(1 - z1) < err1) % z1 == 1
            NUCR = NUCR + 1;
            N_1 = N_1 + 1;
         end
      else % z1 is a non-real unit circle root.
         NUCR = NUCR + 2; % The root and its conjugate.
         ixUCR = (1 + ceil(N*angl(z1)/(2*pi))) : Nspectrum;
         amp(ixUCR) = -amp(ixUCR); % Unwind the amplitude spectrum.
      end
   else
      if(absZ1 < 1)
         inside = 1; % Is this root inside the unit circle?
      else
         inside = 0;
         z1 = 1/z1;
         imagZ1 = imag(z1);
      end
      realZ1 = real(z1);
      
      % The spectrum of [f(n) n=1,2,...,N] is:
      % spectrum(j) = sum(f(k)*exp(-2*pi*k*j/N)) The sum is over k=0,1,...,(N-1)
      % A real root, r, is the linear factor f(n) = [1, -r, 0, 0, 0, ...] so
      % spectrum(j) = f(0)*exp(0) + f(1)*exp(-2*pi*j/N) + 0 + 0 ...
      %             = 1*1 -r*exp(-2*pi*j/N)
      % A complex root, r, is the quadratic factor f(n)=[1,-2*real(r),abs(r)^2,0,0,...]
      % spectrum(j) = f(0)*exp(0) + f(1)*exp(-2*pi*j/N) + f(2)*exp(-4*pi*j/N) + 0 + 0 + ...
      %             = 1*1 -2*real(r)*exp(-2*pi*j/N) + abs(r)^2*exp(-4*pi*j/N)
       if(~imagZ1)
          spectrum = 1 - z1*trigTable;
       else
          spectrum = 1 - trigTable.*(2*realZ1 - trigTable*(realZ1*realZ1+imagZ1*imagZ1));
       end
      
      if(inside)
         ph = ph + angl(spectrum);
      else
         ph = ph - angl(spectrum);
         if(isreal(z1))
            NCxMin = NCxMin + 1;
         else
            NCxMin = NCxMin + 2;
         end
      end
   end
end

% At this point, ph begins and ends at 0. Now add in the linear
% phase shift corresponding to unit circle roots and minphase roots
% and the static phase shift corresponding to the multiplicity of 1.

if(ISODD)
   halfPtShift = pi*(0 : (N+1)/2)/N;
else
   halfPtShift = pi*(0 : N/2)/N;
end

ph = ph - (2*NCxMin + NUCR)*halfPtShift - N_1*pi/2;

if(ISODD)
   lastIx = Nspectrum - 1;
   NyquistPhase = mean(ph((Nspectrum-1):Nspectrum));
else
   lastIx = Nspectrum;
   NyquistPhase = ph(Nspectrum);
end
amp = [amp(1:lastIx), ((-1)^N_Minus_1)*amp((Nspectrum-1):-1:2)];
ph =  [ph(1:lastIx),  2*NyquistPhase-ph((Nspectrum-1):-1:2)];

% ph(1) = -N_1*pi/2. If N_1=2, ph(1)=pi.  If N_1=3, ph(1)=1.5*pi.
% In either case we can subtract pi from every value. Etc.
ph = ph + floor(N_1 / 2)*pi;

y = ifft(amp.*complex(cos(ph),sin(ph)));

max1 = max(abs(x+y));
max2 = max(abs(x-y));
reconError = min(max1,max2); % Error reconstructing polynomial.

% x and -x have the same roots so when we unwind amplitude and phase we
% will not know if we will get x or -x. Check and adjust amp accordingly.
if(max1 < 1E-3)
   amp = -amp;
end

[ph, err] = adjustUnwrapping(ph, amp, wrapped_ph);
amp = reshape(amp,sizeX);
ph = reshape(ph,sizeX);

fprintf('\t\tMax correction-multiple error from wrapped phase     = %g\n', err);
fprintf('\t\tMax error reconstructing the polynomial from spectra = %g\n', reconError);

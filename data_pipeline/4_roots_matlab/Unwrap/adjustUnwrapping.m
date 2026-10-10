% function [unwrapped_ph, err] = adjustUnwrapping(unwrapped_ph, amp, wrapped_ph)
%
% PURPOSE
%   Assume we have unwrapped the phase by some alternative method; e.g. using the roots,
%   unwrapz, or zero padding, unwrapm & unwrapph. The results should exactly equal the
%   wrapped phase modulo pi. Make adjustments so that this is true.
%
% REQUIRED INPUT ARGUMENTS
%   unwrapped_ph - The phase unwrapped by an alternative method.
%   amp - The amplitude spectrum.
%   wrapped_ph - The wrapped phase.
%
% REQUIRED OUTPUT ARGUMENTS
%   unwrapped_ph - The unwrapped phase after adjustments.
%
% OPTIONAL OUTPUT ARGUMENTS
%   err - An indication of the maximum phase alteration made at any frequency so that the 
%         unwrapped phase obtained from the roots agrees with the traditional wrapped 
%         phase mod pi. 0<=err<=1.
%
%         0 - Complete success. Not even one bit of round-off error in any phase value.
%             For 100-degree seismic data, expect values like 1E-13. For 1,000-degree
%             seismic data, expect values like 1E-10.
%         1 - Total failure. Anything close to 1 indicates that the original unwrapping
%             was extremely poor.
% AUTHOR
%   James W. Fox 1/2/04
%
% It is assumed that the phase was unwrapped using some alternative mehtod; e.g. using
% the roots. Assume the phase unwrapping was a "success." You would hope that each value 
% in (unwrapped phase - wrapped phase) would be almost exactly an integral multiple of pi.
%
% I.e., assume at some frequency we have
%    unwrapped phase from roots = wrapped phase + 14.00000000015334*pi
% Obviously the correct value is 14*pi. Therefore, we can correct our unwrapped
% phase at this frequency by replacing it with (wrapped phase + 14*pi). This
% also ensures the wrapped and unwrapped phase agree mod pi. This procedure 
% automatically corrects small unwrapped phase errors caused by slightly 
% incorrect root locations.

% However, there is something to watch out for. If an n'th root of unity
% is a root of the polynomial, the classical wrapped spectrum may be
% incorrect. Why? The FFT is supposed to be exactly 0 at the frequency 
% corresponding to the n'th root of unity. Matlab's angle() will report
% that the phase is 0. However, you can use any random value for the phase at
% that frequency and you will still get the correct answer if you combine the
% amplitude and phase spectra and inverse transform. The correct unwrapped phase
% value is not 0 and it is not a multiple of pi. Whatever the correct value is,
% it is killed by the zero in the amplitude spectrum. Or even worse, it may not
% be killed. A polynomial was created with an n'th root of unity as one of its
% roots. Due to round-off errors in the FFT, the spectrum at that frequency was
% not 0. It was 2.77E-16 - 1.1102E-16 which angle() reported as a phase value of
% -0.38. This was incorrect. Therefore we should not apply the procedure 
% discussed in the previous paragraph if we are at a frequency and there is
% a root of the polynomial that is extrememly close to the corresponding
% n'th root of unity. I.e., do not do this adjustment if the corresponding 
% amplitude spectral value is unusually small.

% I have no evidence that my test for "not being extremely close to an n'th root
% of unity" is an adequate test.

function [unwrapped_ph, err] = adjustUnwrapping(unwrapped_ph, amp, wrapped_ph)

amp = abs(amp); % In case they used an unwrapped amplitude spectrum.
k = find(amp > 1e-13*max(amp)); % Stable frequencies; no nearby root.
multiple = (unwrapped_ph - wrapped_ph)/pi; % Hopefully, virtually an integer array.
roundedMultiple = round(multiple(k));
unwrapped_ph(k) = roundedMultiple*pi + wrapped_ph(k);

err = multiple(k) - roundedMultiple; % Ideally, all values are virtually 0.
% max(abs(err)) is at most 0.5. Multiply by 2 so that 1 indicates total failure.
err = 2*max(abs(err));


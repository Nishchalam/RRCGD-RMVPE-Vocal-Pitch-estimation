% function [coeff, spectrum] = unfactor(z, polyType, Flag, Nspectrum);
%
% PURPOSE:
%   Construct a polynomial from its roots with the highest order coefficient
%   being 1.0. Do the polynomial multiplications in the Fourier domain.
%
% REQUIRED INPUT:
%   z - The roots. 
%       Real roots must truly be real. Eg., if
%       imag(z)==2e-16 it will be considered a non-real root.
%
% OPTIONAL INPUT:
%   polyType - The data type of the output polynomial coefficients.
%          Default = 2.
%          1 - Real coefficients.
%              You only need to supply the upper half-plane roots.
%              Lower half-plane roots will be automatically generated
%              from them. If you supply lower half-plane roots,
%              they will be ignored.
%          2 - Complex coefficients.
%               You must supply all the roots because lower
%               half-plane roots can not be automatically generated
%               from upper half-plane roots.
%
%               unfactor(z) and unfactor(z,2) are equivalent to 
%               Matlab's poly(z). However, they may differ in the least
%               significant digits.
%   Flag - Flag to signal whether to generate coeff.
%          0 - Coeff will be returned as [].
%          1 - Coeff will be non-empty.
%          Default = 1.
%          The algorithm requires that the spectrum be generated but
%          it does not require that coeff be generated.
%          For some applications, we only want the spectrum of
%          the coefficients, but not the coefficients themselves.
%          In this case, setting Flag=0 will save you from an
%          ifft() you do not need.
%   Nspectrum - The desired number of samples in spectrum. Spectrum is
%          probably going to be used to do polynomial division in
%          the Fourier domain. You must ensure that Nspectrum will be
%          the same size as the spectrum of the numerator.
%          Default = 2^nextpow2(1 + the # of roots in both half-planes).
%
% OUTPUT:
%   coeff - Polynomial coefficients. If Flag==0, it will be [].
%   spectrum - The spectrum of coeff.

% USES:
%   lejaJWF, permuteVDC
%
% AUTHOR:
%   James W. Fox 5/1/2001.

% COMMENTS:
% 1) Unfactoring the roots of the Wilkinson polynomial, 1:20, is better carried
% out in the time domain.
% 2) Consider x=[1 0 0 0 1e-100] which has upper half plane roots:
% sqrt(2)/2 + sqrt(2)/2*i and sqrt(2)/2 - sqrt(2)/2*i. roots and lroots do a very
% good job of finding the roots. However, neither unfactor nor poly does a good job.
% 3) unfactor can underflow to give a polynomial where the constant term is 0 even
% though 0 is not a root, (but sqrt2)/2+sqrt(2)/2*i is a root).

function [coeff, spectrum] = unfactor(z, polyType, Flag, Nspectrum);

if(nargin < 1)
  error('Supply $1, the array of roots.');
elseif(nargout > 2)   
  error('Supply at most 2 output arguments.');
elseif(isempty(z))
  coeff = [];
  spectrum = [];
  return;
end

if(nargin < 2)
   polyType = 2;
elseif(length(polyType) ~= 1)
  error('$2 must be a scalar');
elseif(polyType ~= 1 & polyType ~= 2)
  error('$2, must equal 1 or 2.');
end
realPoly = (polyType == 1);

if(nargin < 3)
  Flag = 1;
elseif(length(Flag) ~= 1)
  error('$3 must be a scalar');
elseif(Flag ~= 0 & Flag ~= 1)
   error('$3 must equal 0 or 1.');
end

if(realPoly)
   z = hp2fp(z(find(imag(z) >= 0))); % Ensure it is exactly conjugate symmetric.
end
degree = length(z);

if(nargin < 4)
   Nspectrum = 2^nextpow2(degree + 1);
elseif(Nspectrum ~= 2^nextpow2(Nspectrum))
  error('$4 must be a power of 2.');
elseif(Nspectrum < (degree + 1))
  error('$4 is too small for this many roots.');
end

if(length(z) == 1)
  if(imag(z) == 0 | ~realPoly)
    coeff = [1, -z];
    % Due to recent changes, the following 2 lines will never be accessed. The next 2 lines
    % were designed for the case where upon entry z is a single complex number with non-zero
    % imaginary part and it is desired to create a real polynomial.. However, the 
    % previous call to hp2fp will ensure that its complex conjugate is also added to z. 
    % Thus length(z == 1) will be false.
  else
    absR = abs(z);
    coeff = [1, -2*real(z), absR*absR];
  end
  spectrum = fft(coeff, Nspectrum);
  return;
elseif(isreal(z) | degree < 500)
  coeff = poly(lejaJWF(z));
  spectrum = fft(coeff, Nspectrum);
  return;
end

z = z(:); % Matlab's roots() returns a column vector.
% If you do not do the following, this function can overflow
% for degree 5,000 and higher.
[junk, j] = sort(angle(z));
z = permuteVDC(z(j)); 

% Table of exp(-2*pi*i*j/Nspectrum) for j = 0,1,2,...
if(realPoly)
  m = Nspectrum/2 + 1; % Nyquist index
  k = [0:m-1]*(2*pi/Nspectrum);
else
  k = [0:Nspectrum-1]*(2*pi/Nspectrum);
end
rootsOfUnity = complex(cos(k), -sin(k));
if(realPoly)
   rootsOfUnity(m) = -1; % In Matlab 5.3, cos(pi)=-1 but sin(pi)=-1.2246e-016i
   rootsOfUnity(Nspectrum/4+1) = -i; % In Matlab 5.3, cos(pi/2)=6.1232e-017 but sin(pi/2)=1
end

spectrum = rootSpectrum(z, rootsOfUnity); % The spectrum of the polynomial determined by z.
if(realPoly) % Supply the conjugate half we did not compute.
   spectrum = [spectrum, conj(spectrum(m-1:-1:2))];
   % Due to roundoff errors, things that are supposed to be real may not be.
   spectrum(1) = real(spectrum(1));
   spectrum(m) = real(spectrum(m)); 
end

if(Flag == 0)
  coeff = [];
elseif(Flag == 1)
  coeff = ifft(spectrum);
  coeff = coeff(1 : (degree+1));

  if(realPoly)
     coeff = real(coeff);
  end

  % On one occasion, ifft() above underflowed to zero while computing
  % the highest order coefficient. Thus, this routine returned
  % a n-1 degree polynomial given n roots. In that instance the
  % coefficients were very flat (and small). If this ever happens
  % again, I reset coeff(1) to the first non-zero coefficient.
  if(~coeff(1))
    coeff(1) = coeff(ixFirstNZ(coeff));
  end

  % Matlab's routine poly(), returns monic polynomials.
  coeff = coeff ./ coeff(1);
else
  error('$3 must equal 0 or 1');
end

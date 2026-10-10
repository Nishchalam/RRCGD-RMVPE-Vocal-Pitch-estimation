% function spectrum = rootSpectrum(z, rootsOfUnity)
%
% PURPOSE
%   Given an array of roots, z, and the roots of unity appropriate for this
%   length, find the Fourier spectrum of the polynomial 
%   determined by unfactoring z. This is done by multiplying together the zero
%   padded spectrum associated with each root.
%
% REQUIRED INPUT:
%   z = The roots
%   rootsOfUnity - Round up the degree of the polynomial determined by z to a
%                  power of 2. This is the roots of unity for this length.
%
% OUTPUT:
%   spectrum - The spectrum of the polynomial determined by z.

% USES:
%   nothing
%
% AUTHOR:
%   James W. Fox 1/1/2006

% The spectrum of f(n) n=1,2,... is:
% spectrum(j) = sum(f(k)*exp(-2*pi*k*j/N)) The sum is over k=0,1,...
% The linear factor z-r corresponds to f(n) = [1, -r, 0, 0, 0, ...] so
% spectrum(j) = f(0)*exp(0) + f(1)*exp(-2*pi*j/N) + 0 + 0 ...
%             = 1*1 -r*exp(-2*pi*j/N)
% This works even for complex r.

function spectrum = rootSpectrum(z, rootsOfUnity)

global DIAGNOSTICS nCPUs

if(isempty(DIAGNOSTICS))
  DIAGNOSTICS = 0;
  nCPUs = 1;  
end

spectrum = 1;
if(nCPUs == 1) % There is only 1 cpu.
  for k = 1:length(z)
    if(DIAGNOSTICS) % Useful for really high degree polynomials to let you know how much is left to do.
      if(~rem(k,1000))
        fprintf(1,'\t Unfactoring root # %d\n', k);
      end
    end
    spectrum = spectrum .* (1 - z(k)*rootsOfUnity); % Evaluate 1-z(k)*z at the roots of unity. This is the FFT.
  end
else
   % Partition z into the first 1/nCPUs roots, and the second 1/nCPUs roots etc. Then send them to the nCPUs different
   % CPUs and do the above computation nCPUs times to get multiple partial spectra. Then multiply together the spectra.
   % See bigPolish for my suggestion on how to do this.
   %
   % Note, it is extremely important that you partition z into the first 1/nCPUs values and not z(1:nCPUs:end). 
   % Doing this could cause disasterous results, depending on the length.
end


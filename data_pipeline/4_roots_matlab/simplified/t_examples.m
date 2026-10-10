% Test minimum phase conversion, etc.
%
function t_examples(n);
% Default n?
if (nargin < 1)
    n = 2001;
 end
 
fprintf('\n\nNine plots will be drawn. Please be patient.\n');
fprintf('The second one takes a very long time to create.\n');
 
% Initialize.
dt = 0.004; fl = 5; fu = 25; 
fn = 1/(2*dt); bw = 1 - (fu-fl)/fn;
tstr = ['Time in Secs. (dt = ', num2str(round(dt*1000)), ' MSec.)'];
load s.dat
figure(1); plot([0:n-1]*dt, s);
axis([0, (n-1)*dt, -1.2, 1.2]);
title(sprintf('Linear Chirp: %d Hz. to %d Hz.', fl, fu));
xlabel(tstr);
ylabel('Chirp Amplitude');
drawnow;

% Minimum phase.
[ms, zo, zm] = rminphw(s);
nn = round(n*bw/5); tms = ms(1:nn);
figure(2); plot([0:nn-1]*dt, tms);
axis([0, (nn-1)*dt, min(tms), max(tms)]);
title('Minimum Phase Equivalent of Chirp');
xlabel(tstr);
ylabel('Response Amplitude');
disp(sprintf('Spectrum error = %e', ...
   norm(abs(fft(s))-abs(fft(ms)), inf)));
drawnow;

% Unwrap phase.
[ps, err] = runwrap(s, 1e-12, zo);
disp(sprintf('Phase error = %e', err(1)));
m = fix(length(ps)); df = 1/((m-1)*2*dt);
tps = ps(1:m);
figure(3); plot([0:m-1]*df, tps);
axis([0, (m-1)*df, min(tps), max(tps)]);
title(sprintf('Chirp, Root UW-Phase: %d Hz. to %d Hz.', fl, fu));
xlabel('Frequency in Hz.');
ylabel('UW-Phase in Radians');
drawnow;

% MP Root unwrap phase.
pms = runwrap(ms, 1e-8, zm);
tpms = pms(1:m);
figure(4); plot([0:m-1]*df, tpms);
axis([0, (m-1)*df, min(tpms), max(tpms)]);
title(sprintf('Chirp (Min. Ph. Eq.), Root UW-Phase: %d Hz. to %d Hz.', fl, fu));
xlabel('Frequency in Hz.');
ylabel('UW-Phase in Radians');
drawnow;

% MP Fourier unwrap phase.
pms = unwrap(angle(fft(ms))); 
tpms = pms(1:m);
figure(5); plot([0:m-1]*df, tpms);
axis([0, (m-1)*df, min(tpms), max(tpms)]);
title(sprintf('Chirp (Min. Ph. Eq.), Fourier UW-Phase: %d Hz. to %d Hz.', fl, fu));
xlabel('Frequency in Hz.');
ylabel('UW-Phase in Radians');
drawnow;

% Fourier unwrap phase.
pms = unwrap(angle(fft(s)))'; 
tpms = pms(1:m);
figure(6); plot([0:m-1]*df, tpms);
axis([0, (m-1)*df, min(tpms), max(tpms)]);
title(sprintf('Chirp, Fourier UW-Phase: %d Hz. to %d Hz.', fl, fu));
xlabel('Frequency in Hz.');
ylabel('UW-Phase in Radians');
drawnow;

% Sweep spectrum.
fs = abs(fft(s))';
tfs = fs(1:m);
figure(7); plot([0:m-1]*df, tfs);
axis([0, (m-1)*df, min(tfs), max(tfs)]);
title(sprintf('Chirp Amplitude Spectrum: %d Hz. to %d Hz.', fl, fu));
xlabel('Frequency in Hz.');
ylabel('Fourier Magnitude');
drawnow;

% Roots.
figure(8); proots(1.0./zo);

figure(9); proots(1.0./zm);
%
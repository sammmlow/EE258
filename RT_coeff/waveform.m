%% GP 258J: Simulated waveform

% Free space EM properties
eps0 = 8.854e-12;
mu0 = 1.257e-6;
c0 = 1 / sqrt(eps0*mu0);
Z0 = mu0 * c0;

% Cole-Cole model parameters
eps_inf = 3.140;
eps_delta = 1.708;
tau = 14.65e-12;
alpha = 0.061;
sigma_s = 0.036;

% Decibel level
dBp = @(x) 10 .* log10(x);
dBa = @(x) 20 .* log10(x);

%% Waveform distortion

%%% Complex chirp signal %%%

f0 = 12.5e9;    % Carrier (Hz)
T = 10e-9;      % Pulse length (s)
fs = 60e9;      % Sampling freq. (Hz)
B = 10e9;       % Bandwidth (Hz)
achirp = B/T;   % Chirp (FM) modulation

dt = 1 / fs;                    % Time step (s)
t = (0:dt:T-dt) - T/2;          % Time samples (s)
s = exp(1j*2*pi*f0.*t) .* ...
    exp(1j.*pi*achirp.*t.^2);   % Chirp signal
s = s .* tukeywin(length(s), 1/10)';
s_tmp = s;                      % Matched filter

% Padding with zeros
pad_time1 = 5e-9;  pad_time2 = 10e-9;
pad_samp1 = round(pad_time1 * fs);
pad_samp2 = round(pad_time2 * fs);

% Incidence waveform
s = [zeros(1, pad_samp1), s, zeros(1, pad_samp2)];
Nt = length(s);  t = (0:Nt-1) .* dt;

% Chirp spectrum
s_hat = fft(s);

%%% Transfer function %%%

% Frequency samples (Hz)
f = (-Nt/2:Nt/2-1) * (fs / Nt);  f = fftshift(f);
w = 2*pi .*f;

% Wavelength (m)
lambda = c0 ./ f;

% Complex dielectric constant
epsd = eps_inf + eps_delta ./ (1 + (1j*tau.*w).^(1-alpha)) ...
    + sigma_s ./ (1j*eps0.*w);

% Thickness of tissue (m)
h = [1e-2, 2e-2, 5e-2]';

% Refractive index
nd = sqrt(epsd);  kappa = -imag(nd);

% Exponent term
exp_gamma = exp(-1j*2*pi.*h./lambda.*nd);

% Transmission
Ra = (1 - nd) ./ (1 + nd);
HT = (1 - Ra.^2) .* exp_gamma ./ (1 - (Ra.*exp_gamma).^2);

% Reflection
HR = Ra .* (1 - exp_gamma.^2) ./ (1 - (Ra.*exp_gamma).^2);

% Fix zero component
HT(:, 1) = 0;  HR(:, 1) = 0;

%%% Distorted waveform %%%

% Inverse FFT
s2_hat = s_hat .* HT;
s2 = ifft(s2_hat, [], 2);

% Plot signal
figure('Name', 'Waveforms');
ax1 = subplot(4,1,1);  plot(t.*1e9, real(s), 'k-');
title('Transmitted chirp pulse');
ax2 = subplot(4,1,2);  plot(t.*1e9, real(s2(1, :)), 'k-');
title(sprintf('h = %.1f cm', h(1).*1e2));
ax3 = subplot(4,1,3);  plot(t.*1e9, real(s2(2, :)), 'k-');
title(sprintf('h = %.1f cm', h(2).*1e2));
ax4 = subplot(4,1,4);  plot(t.*1e9, real(s2(3, :)), 'k-');
title(sprintf('h = %.1f cm', h(3).*1e2));
linkaxes([ax1, ax2, ax3, ax4], 'x');
xlabel('Time (ns)');  % xlim([5, 16]);

%% Matched filter detection

% Add noise
noise_power = 3.0 .* std(s2, [], 2);
s2 = s2 + sqrt(noise_power./2) .* (randn(1, Nt) + 1j*randn(1, Nt));

% Scaling factor
A = sum(abs(s_tmp).^2);     

% Matched filter
s_filt = conv(s, conj(fliplr(s_tmp)), 'same') ./ A;
s2_filt = zeros(size(s2));
for i = 1:length(h)
    s2_filt(i, :) = conv(s2(i, :), conj(fliplr(s_tmp)), 'same') ./ A;
end

% Plot signal
figure('Name', 'Waveforms');
ax1 = subplot(4,1,1);  plot(t.*1e9, real(s_filt), 'k-');
title('Transmitted chirp pulse');
ax2 = subplot(4,1,2);  plot(t.*1e9, real(s2_filt(1, :)), 'k-');
title(sprintf('h = %.1f cm', h(1).*1e2));
ax3 = subplot(4,1,3);  plot(t.*1e9, real(s2_filt(2, :)), 'k-');
title(sprintf('h = %.1f cm', h(2).*1e2));
ax4 = subplot(4,1,4);  plot(t.*1e9, real(s2_filt(3, :)), 'k-');
title(sprintf('h = %.1f cm', h(3).*1e2));
linkaxes([ax1, ax2, ax3, ax4], 'x');
xlabel('Time (ns)');  % xlim([5, 16]);

%% GP 258J: Simulated waveform

%%% Free space EM properties %%%

eps0 = 8.854e-12;           % Permittivity (F/m)
mu0  = 1.257e-6;            % Permeability (H/m)
c0   = 1 / sqrt(eps0*mu0);  % Speed of light (m/s)
Z0   = mu0 * c0;            % Impedance (Ω)


%%% Cole-Cole model parameters %%%

materials = struct();

% Low water content: > 85 % adipose
materials.lwc = struct(...
    'eps_inf',   3.140, ...     % High-frequency permittivity
    'eps_delta', 1.708, ...     % Permittivity dispersion
    'tau',       14.65e-12, ... % Relaxation time (s)
    'alpha',     0.061, ...     % Distribution parameter
    'sigma',     0.036 ...      % Static conductivity (S/m)
);

% Medium water content: 30 % - 85 % adipose
materials.mwc = struct(...
    'eps_inf',   5.573, ...
    'eps_delta', 34.57, ...
    'tau',       9.149e-12, ...
    'alpha',     0.095, ...
    'sigma',     0.524 ...
);

% High water content: < 30 % adipose
materials.hwc = struct(...
    'eps_inf',   7.821, ...
    'eps_delta', 41.48, ...
    'tau',       10.66e-12, ...
    'alpha',     0.047, ...
    'sigma',     0.713 ...
);

% Add free space properties
materialNames = fieldnames(materials);
for i = 1:numel(materialNames)
    matName = materialNames{i};
    materials.(matName).eps0 = eps0;
    materials.(matName).mu0  = mu0;
    materials.(matName).c0   = c0;
    materials.(matName).Z0   = Z0;
end


%%% Conversion to decibel level %%%

dBp = @(x) 10 .* log10(x);
dBa = @(x) 20 .* log10(x);


%% Waveform distortion

% Select material type
mat_sel = 'mwc';


%%% Complex chirp signal %%%

f0 = 12.5e9;    % Carrier (Hz)
T = 10e-9;      % Pulse length (s)
fs = 64e9;      % Sampling freq. (Hz)
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
Nt = length(s);  t = (0:Nt-1) .* dt - pad_time1;

% Chirp spectrum
s_hat = fft(s);


%%% Transfer function %%%

% Frequency samples (Hz)
f = (-Nt/2:Nt/2-1) * (fs / Nt);  f = fftshift(f);
w = 2*pi .*f;

% Wavelength (m)
lambda = c0 ./ f;

% Complex dielectric constant
for i = 1:numel(materialNames)
    matName = materialNames{i};
    epsd.(matName) = cole_model(w, materials.(matName));
end

% Thickness of tissue (m)
h = [1e-2, 2e-2, 5e-2]';

for i = 1:numel(materialNames)
    matName = materialNames{i};
    
    nd_ = sqrt(epsd.(matName));     % Complex refractive index
    R_ = (1 - nd_) ./ (1 + nd_);    % Reflection coeff. (air -> tissue)

    % Transfer function
    lambda = c0 ./ f;  
    exp_gamma = exp(-1j*2*pi.*h./lambda.*nd_);
    HT_ = (1-R_.^2) .* exp_gamma ./ (1 - (R_.*exp_gamma).^2);
    HR_ = R_ .* (1 - exp_gamma.^2) ./ (1 - (R_.*exp_gamma).^2);

    % Fix zero component
    HT_(:, 1) = 0;  HR_(:, 1) = 0;

    % Store results
    nd.(matName) = nd_;  R.(matName) = R_;
    HT.(matName) = HT_;  HR.(matName) = HR_;
    Tp.(matName) = Tp_;  Rp.(matName) = Rp_;

end


%%% Distorted waveform %%%

% Inverse FFT
s2_hat = s_hat .* HT.(mat_sel);
s2 = ifft(s2_hat, [], 2);

% Plot signal
figure('Name', 'Waveforms');
ax1 = subplot(4,1,1);  plot(t.*1e9, real(s), 'k-');  hold on;
plot(t.*1e9, abs(s), 'r-');
title('Emitted chirp pulse');
ax2 = subplot(4,1,2);  plot(t.*1e9, real(s2(1, :)), 'k-');  hold on;
plot(t.*1e9, abs(s2(1, :)), 'r-');
title(sprintf('Transmission, h = %.1f cm', h(1).*1e2));
ax3 = subplot(4,1,3);  plot(t.*1e9, real(s2(2, :)), 'k-');  hold on;
plot(t.*1e9, abs(s2(2, :)), 'r-');
title(sprintf('Transmission, h = %.1f cm', h(2).*1e2));
ax4 = subplot(4,1,4);  plot(t.*1e9, real(s2(3, :)), 'k-');  hold on;
plot(t.*1e9, abs(s2(3, :)), 'r-');
title(sprintf('Transmission, h = %.1f cm', h(3).*1e2));
linkaxes([ax1, ax2, ax3, ax4], 'x');
xlabel('Time (ns)');  xlim([-0.5, 10.5]);


%% Matched filter detection

% Add noise
for i = 1:length(h)
    s2(i,:) = awgn(s2(i,:), 1);
end

% Min refractive index
n_ref = min(real(nd.(mat_sel)(2:end)));

% Scaling factor
A = sum(abs(s_tmp).^2);     

% Matched filter
t_lag = t(1:Nt-T*fs+1);
s_filt = conv(s, conj(fliplr(s_tmp)), 'valid') ./ A;
s2_filt = zeros(size(t_lag));
for i = 1:length(h)
    s2_filt(i, :) = conv(s2(i, :), conj(fliplr(s_tmp)), 'valid') ./ A;
end

% Plot signal
figure('Name', 'Waveforms');
ax1 = subplot(4,1,1);  plot(t_lag.*1e9, abs(s_filt), 'k-');
xline(0, 'r-', 'LineWidth', 2);
title('Emitted chirp pulse');
ax2 = subplot(4,1,2);  plot(t_lag.*1e9, abs(s2_filt(1, :)), 'k-');
xline(h(1)/c0*n_ref*1e9, 'r-', 'LineWidth', 2);
title(sprintf('Transmission, h = %.1f cm', h(1).*1e2));
ax3 = subplot(4,1,3);  plot(t_lag.*1e9, abs(s2_filt(2, :)), 'k-');
xline(h(2)/c0*n_ref*1e9, 'r-', 'LineWidth', 2);
title(sprintf('Transmission, h = %.1f cm', h(2).*1e2));
ax4 = subplot(4,1,4);  plot(t_lag.*1e9, abs(s2_filt(3, :)), 'k-');
% xline(h(3)/c0*n_ref*1e9, 'r-', 'LineWidth', 2);
title(sprintf('Transmission, h = %.1f cm', h(3).*1e2));
linkaxes([ax1, ax2, ax3, ax4], 'x');
xlabel('Time (ns)');  xlim([-0.25, 0.75]);

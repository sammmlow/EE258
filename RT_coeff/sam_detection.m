%% Copy of Jiqing's code from "waveform.m"

clc; clear all; close all;
set(groot, 'DefaultLineLineWidth', 1);
defStrFig = 'defaultFigureUnits';
defStrUnits = 'inches';
defStrPos = 'defaultFigurePosition';
set(groot, defStrFig, defStrUnits, defStrPos, [0.5 0.5 7.0 7.5]);

%% Set the SNR to loop through.
SNR_ALL = [100];

% Set the transmitter-to-lesion and lesion-to-receiver distances.
Rt = 0.05;
Rr = 0.05;

% Free space EM properties
eps0 = 8.854e-12;
mu0 = 1.257e-6;
c0 = 1 / sqrt(eps0*mu0);

f0 = 2.4e9;   % Carrier (Hz)
T = 1e-9;      % Pulse length (s)
fs = 100e9;    % Sampling freq. (Hz)
B = 10e9;     % Bandwidth (Hz)
achirp = B/T;  % Chirp (FM) modulation

h = 0.05; % Thickness of tissue (m)

pad_time = 4e-9; % Duration of zero pad after chirp

% Set time step and the initial time axis for samples within pulse.
dt = 1 / fs;       
ti = (0:dt:T-dt);  % Initial time axis for chirp samples (s)

% Generate noiseless chirp pulse and matched filter
s = exp(1j*2*pi*f0.*ti) .* exp(1j.*pi*achirp.*ti.^2);
s_match = conj( fliplr( s ) ); % Matched filter

% Add the `tukeywin` after generating the matched filter.
s = s .* tukeywin(length(s), 1/10)';

% Pad signal with zeros
pad_samples = round(pad_time * fs);

% Incidence waveform
s = [ s, zeros(1, pad_samples) ];
Nt = length(s);

% Frequency axis
f = ( -Nt/2 : Nt/2-1 ) * (fs / Nt);
f = fftshift(f);
w = 2*pi .*f;
w0 = 2*pi * f0;

% Time axis
t = (0:Nt-1) .* dt;

% Wavelength (m)
lambda = c0 ./ f;

% Cole-Cole model parameters
eps_inf_all   = [3.140, 5.573, 7.821];
eps_delta_all = [1.708, 34.57, 41.48];
tau_all       = [14.65e-12, 9.149e-12, 10.66e-12];
alpha_all     = [0.061, 0.095, 0.047];
sigma_s_all   = [0.036, 0.524, 0.713];

model = ["Low water: ", ...
         "Medium water: ", ...
         "High water: "];

colors = {[0.4350 0.0780 0.0840], ... % Red 1
          [0.9350 0.3780 0.3240], ... % Red 2
          [0.0000 0.2470 0.5410], ... % Blue 1
          [0.3000 0.4470 0.9810]};    % Blue 2

% Loop through each model parameter for low and medium water content.
for k = 1:2

    % Set the color according to the tissue type.
    if k == 1
        colorIdx = 1;
    elseif k == 2
        colorIdx = 3;
    end
    
    % Tissue model parameters
    eps_inf = eps_inf_all(k);
    eps_delta = eps_delta_all(k);
    tau = tau_all(k);
    alpha = alpha_all(k);
    sigma_s = sigma_s_all(k);

    model_str = char(model(k));

    % Lesion model parameters
    eps_inf_tumor = eps_inf_all(3);
    eps_delta_tumor = eps_delta_all(3);
    tau_tumor = tau_all(3);
    alpha_tumor = alpha_all(3);
    sigma_s_tumor = sigma_s_all(3);

    % Complex dielectric constant (across all frequencies)
    epsd = eps_inf + eps_delta ./ (1 + (1j*tau.*w).^(1-alpha)) ...
        + sigma_s ./ (1j*eps0.*w);
    
    % Complex dielectric constant at the carrier frequency,
    epsd_f0 = eps_inf + eps_delta ./ (1 + (1j*tau*w0).^(1-alpha)) ...
        + sigma_s ./ (1j*eps0*w0);

    % Complex dielectric constant of the tumour at the carrier frequency,
    epsd_f0_tumour = eps_inf_tumor + eps_delta_tumor ./ ...
        (1 + (1j*tau_tumor*w0).^(1-alpha_tumor)) ...
        + sigma_s_tumor ./ (1j*eps0*w0);
    
    % Refractive index
    nd = sqrt(epsd);
    nd0 = sqrt(epsd_f0);
    kappa = -imag(nd);  % UNUSED

    % Compute the propagation constant.
    gamma = 1j * w0 * sqrt(mu0 * eps0 * epsd_f0);
    attenuation = exp(-( real(gamma) * (Rr + Rt) ));
    reflection_coeff = (sqrt(epsd_f0_tumour) - sqrt(epsd_f0)) / ...
        (sqrt(epsd_f0_tumour) + sqrt(epsd_f0));
    
    % Exponent term in the transfer function (for main signal)
    exp_gamma = exp(-1j*2*pi.*h./lambda.*nd);
    
    % Transmission transfer function (includes internal reflections)
    Ra = (1 - nd) ./ (1 + nd);
    HT = (1 - Ra.^2) .* exp_gamma ./ (1 - (Ra.*exp_gamma).^2);
    HT(1) = 0; % Fix zero component
    
    % For each SNR, try to generate the corresponding waveform 
    for i = 1 : length(SNR_ALL)

        id = 10 * i + 100 * k; % Convenient number for identifying plots
        idStr = num2str(id);
    
        SNR = SNR_ALL(i);
    
        % Add AWGN noise to the chirp signal
        s_noisy = awgn( s, SNR);
        
        % Chirp spectrum
        s_fft = fft(s_noisy);
        s_fft_output = s_fft .* HT;
        s_received = ifft(s_fft_output, [], 2);

        % Model the reflected signal from the lesion
        s_reflected = reflection_coeff * attenuation * s_noisy;

        % Model the delay time of the reflected signal by taking the
        % Fourier Transform of the original signal, then apply a linear
        % phase shift (corresponding to time delay). Then take the Inverse
        % Fourier Transform.
        delay_time = (Rt + Rr) / (c0 / real(nd0));
        delay_samples = round(delay_time * fs);
        s_reflected_delayed = [ zeros(1, delay_samples) ...
            s_reflected(1:(end-delay_samples)) ];

        % Add the reflected signal to the received signal
        s_received_plus_reflected = s_received + s_reflected_delayed;
        
        % Plot signal of transmitted chirp
        fig1 = figure(id + 1);
        subplot(2,1,1);  plot(t.*1e9, real(s_noisy), 'k', ...
            'DisplayName', sprintf([model_str ...
            'Transmitted chirp pulse at SNR = %.1f dB'], SNR)); 
        title('Transmitted chirp pulse');
        xlabel('Time (ns)'); legend show; hold on; grid on;
        legend('Location', 'southeast');
        
        % Plot the received signal against a distance axis (x)
        subplot(2,1,2);
        plot(t.*1e9, real(s_received), ...
            'DisplayName', sprintf([model_str ...
            'Received (no tumour), SNR = %.1f dB'], SNR), ...
            Color=colors{colorIdx});
        hold on; grid on;
        plot(t.*1e9, real(s_received_plus_reflected), ...
            'DisplayName', sprintf([model_str ...
            'Received (with tumour), SNR = %.1f dB'], SNR), ...
            Color=colors{colorIdx+1});
        title(sprintf(['Received signal: h = %.1f cm, SNR = %.1f dB'], ...
            h.*1e2, SNR));
        xlabel('Time (ns)'); legend show;
        legend('Location', 'southeast');
    
        % Export graphics.
        exportgraphics(fig1, ...
            ['../Figures/detection_signal' idStr ...
            'SNR=' num2str(SNR) '.png'], ...
            'Resolution', 300);
        
        % To compute distance travelled. Light slows down under refraction.
        % dist = t .* c0 / real(nd0);  % UNUSED: 
        
        % Matched filter output for the original transmitted signal
        s_filtered = conv(s, s_match, 'same');

        % Matched filter output for the noisy output signal (NO tumour)
        s_filtered_noisy = conv(s_received, s_match, 'same');
        
        % Matched filter output for the noisy output signal (with tumour)
        s_filtered_noisy_plus_reflected = ...
            conv(s_received_plus_reflected, s_match, 'same');
    
        % Plot match filtered output of the transmitted chirp
        fig2 = figure(id + 2);
        subplot(2,1,1);
        plot(t.*1e9, real(s_filtered), 'k-'); hold on; grid on;
        title('Match filter output of transmitted chirp pulse');
        xlabel('Time (ns)'); ylabel('Matched filter output');

        % Plot match filtered output of the noisy received signal
        subplot(2,1,2); 
        plot(t.*1e9, real(s_filtered_noisy), ...
            'DisplayName', sprintf([model_str ...
            'Match filter output (no tumour), SNR = %.1f dB'], SNR), ...
            Color=colors{colorIdx});
        hold on; grid on;
        plot(t.*1e9, real(s_filtered_noisy_plus_reflected), ...
            'LineStyle', '--', 'DisplayName', sprintf([model_str ...
            'Match filter output (with tumour), SNR = %.1f dB'], SNR), ...
            Color=colors{colorIdx+1});
        title(sprintf(['Match filter output of received ' ...
            'signal, h = %.1f cm, SNR = %.1f dB'], h(1).*1e2, SNR));
        xlabel('Time (ns)'); ylabel('Matched filter output'); legend show;
        legend('Location', 'southeast');

        % % Compute what is the expected time of flight of the direct LOS
        % % signal and the reflected signal from the tumour.
        % tof_direct = T/2 + h .*1e9 / (c0 / real(nd0));
        % tof_reflect = T/2 + (Rt + Rr) .*1e9 / (c0 / real(nd0));
        % xline(tof_direct, '--k', 'Expected ToA (direct)', ...
        %     'HandleVisibility', 'off');
        % xline(tof_reflect, '--k', 'Expected ToA (lesion)', ...
        %     'HandleVisibility', 'off');

        % Export graphics.
        exportgraphics(fig2, ...
            ['../Figures/detection_match_filter' idStr ...
            'SNR=' num2str(SNR) '.png'], 'Resolution', 300);
    
    end

end

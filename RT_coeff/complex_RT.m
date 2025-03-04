%% GP 258J: Complex dielectric material

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

%% Dielectric property

% Frequency-dependent property
f = linspace(0.5, 20, 100) * 1e9;  w = 2*pi .* f;
epsd = eps_inf + eps_delta ./ (1 + (1j*tau.*w).^(1-alpha)) ...
    + sigma_s ./ (1j*eps0.*w);

% Wavelength
lambda = c0 ./ f;

figure('Name', 'Dielectric');
ax1 = subplot(1,3,1);  plot(f./1e9, real(epsd), 'k-');
xlabel('Frequency (GHz)');  title('Rel. Dielectric Constant');
ylabel('$\epsilon^*_{R}$', 'Interpreter', 'latex', 'FontSize', 20);
ax2 = subplot(1,3,2);  plot(f./1e9, -imag(epsd), 'k-');
xlabel('Frequency (GHz)');  title('Rel. Dielectric Loss');
ylabel('$\epsilon^*_{I}$', 'Interpreter', 'latex', 'FontSize', 20);
linkaxes([ax1, ax2], 'x');
ax3 = subplot(1,3,3);  plot(real(epsd), -imag(epsd), 'k-');
xlabel('$\epsilon^*_{R}$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\epsilon^*_{I}$', 'Interpreter', 'latex', 'FontSize', 20);

% ax2 = subplot(1,2,2);  plot(f./1e9, -w.*eps0.*imag(epsd), 'k-');
% xlabel('Frequency (GHz)');  title('Effective Conductivity');
% ylabel('$\sigma$ (S/m)', 'Interpreter', 'latex', 'FontSize', 20);

%% Amplitude RT coeff. (Normal incidence)

% Thickness of tissue (m)
h = [1e-2, 1e-1, 5e-1]';

% Refractive index
nd = sqrt(epsd);  kappa = -imag(nd);

% RT coeff. (Normal incidence: Free space -> Tissue)
Ra = (1 - nd) ./ (1 + nd);
Ta = Ra + 1;

% Absorption (optical length)
A0 = exp(-2*pi*h./lambda.*kappa);

figure('Name', 'RT coeff. (Amp.)');
ax1 = subplot(1,3,1);  plot(f./1e9, abs(Ra), 'k-');
xlabel('Frequency (GHz)');  title('Reflection');
ax2 = subplot(1,3,2);  plot(f./1e9, abs(Ta), 'k-');
xlabel('Frequency (GHz)');  title('Transmission');
ax3 = subplot(1,3,3);  plot(f./1e9, A0, 'k-');
xlabel('Frequency (GHz)');  title('Absorption');
linkaxes([ax1, ax2, ax3], 'x');

%% Apparent power RT coeff. (Normal incidence)

% Power RT coeff.
Ap0 = A0.^2;  Rp0 = abs(Ra).^2;  Tp0 = abs(Ta).^2;

% RT coeff. from infinite sum
Tp_all = Ap0 .* (1 - Rp0).^2 ./ (1 - (Rp0.*Ap0).^2);
Rp_all = (1 + Tp_all.*Ap0) .* Rp0;

figure('Name', 'Apparent RT coeff. (Power)');
ax1 = subplot(1,3,1);  plot(f./1e9, dBp(Rp_all));
xlabel('Frequency (GHz)');  title('Reflection');
ax2 = subplot(1,3,2);  plot(f./1e9, dBp(Tp_all));
xlabel('Frequency (GHz)');  title('Transmission');
ax3 = subplot(1,3,3);  plot(f./1e9, dBp(1-Rp_all-Tp_all));
xlabel('Frequency (GHz)');  title('Absorption');
linkaxes([ax1, ax2, ax3], 'x');

legend_arr = {};
for i = 1:length(h)
    legend_arr{i} = sprintf('h = %.1f cm', h(i).*1e2);
end
legend(legend_arr);

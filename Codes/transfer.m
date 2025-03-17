%% GP 258J: Complex dielectric material

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


%% Dielectric property

%%% Compute relative permissivity for each material %%%

f = linspace(0.5, 20, 100) * 1e9;  w = 2*pi .* f;
for i = 1:numel(materialNames)
    matName = materialNames{i};
    epsd.(matName) = cole_model(w, materials.(matName));
end


%%% Plot relative permissivity %%%

fig = figure('Name', 'Dielectric', ...
    'Position', [100 100 900 350], 'Units', 'pixels');
colors = ["r", "k", "b"];

ax1 = subplot(1,2,1);  title(ax1, 'Rel. Dielectric Constant');
xlabel(ax1, 'Frequency (GHz)');
ylabel(ax1, '$\epsilon^*_{R}$', 'Interpreter', 'latex', 'FontSize', 20);
hold(ax1, 'on');
ax2 = subplot(1,2,2);  title(ax2, 'Rel. Dielectric Loss');
xlabel(ax2, 'Frequency (GHz)');
ylabel(ax2, '$\epsilon^*_{I}$', 'Interpreter', 'latex', 'FontSize', 20);
hold(ax2, 'on');

for i = 1:numel(materialNames)
    matName = materialNames{i};
    plot(ax1, f./1e9, real(epsd.(matName)), 'Color', colors(i));
    plot(ax2, f./1e9, -imag(epsd.(matName)), 'Color', colors(i));
end
linkaxes([ax1, ax2], 'x');

legend(ax1, {"Low water content", "Medium water content", "High water content"});
% exportgraphics(fig, '../Figures/cole.png', 'Resolution', 300);


%% Transfer function (Normal incidence)

%%% Calculate transfer function (Fabry-Perot interferometer) %%%

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

    % Transmission with propagation phase removed (ref: fmax)
    HT2_ = HT_ ./ exp(-1j*2*pi*h./lambda.*real(nd_(end)));

    % Power coefficients
    Rp_ = abs(HR_).^2;  Tp_ = abs(HT_).^2;

    % Store results
    nd.(matName) = nd_;  R.(matName) = R_;
    HT.(matName) = HT_;  HR.(matName) = HR_;  HT2.(matName) = HT2_;
    Tp.(matName) = Tp_;  Rp.(matName) = Rp_;

end


%%% Plot transfer function (amplitude & phase) %%%

% Plot thickness
ind_h = 1;
title_str = sprintf('Tissue thickness %d cm', h(ind_h)*1e2);

fig = figure('Name', 'Transmission', ...
    'Position', [100 100 900 350], 'Units', 'pixels');

ax1 = subplot(1,2,1);  title(ax1, title_str);
xlabel(ax1, 'Frequency (GHz)');  ylabel(ax1, 'Transmission (dB)');
hold(ax1, 'on');
ax2 = subplot(1,2,2);  title(ax2, title_str);
xlabel(ax2, 'Frequency (GHz)');  ylabel(ax2, 'Phase (deg)');
hold(ax2, 'on');

for i = 1:numel(materialNames)
    matName = materialNames{i};
    plot(ax1, f./1e9, dBa(abs(HT.(matName)(ind_h,:))), 'Color', colors(i));
    plot(ax2, f./1e9, rad2deg(angle(HT2.(matName)(ind_h,:))), 'Color', colors(i));
end
linkaxes([ax1, ax2], 'x');  ylim(ax1, [-30, 0]);  ylim(ax2, [-150, 150]);

legend(ax1, {"Low water content", "Medium water content", "High water content"}, ...
     'Location', 'best');
% exportgraphics(fig, '../Figures/transfer_T.png', 'Resolution', 300);


fig = figure('Name', 'Reflection', ...
    'Position', [100 100 900 350], 'Units', 'pixels');

ax1 = subplot(1,2,1);  title(ax1, title_str);
xlabel(ax1, 'Frequency (GHz)');  ylabel(ax1, 'Reflection (dB)');
hold(ax1, 'on');
ax2 = subplot(1,2,2);  title(ax2, title_str);
xlabel(ax2, 'Frequency (GHz)');  ylabel(ax2, 'Phase (deg)');
hold(ax2, 'on');

for i = 1:numel(materialNames)
    matName = materialNames{i};
    plot(ax1, f./1e9, dBa(abs(HR.(matName)(ind_h,:))), 'Color', colors(i));
    plot(ax2, f./1e9, wrapTo360(rad2deg(angle(HR.(matName)(ind_h,:)))), 'Color', colors(i));
end
linkaxes([ax1, ax2], 'x');  ylim(ax2, [120, 240]);

legend(ax1, {"Low water content", "Medium water content", "High water content"}, ...
     'Location', 'best');
% exportgraphics(fig, '../Figures/transfer_R.png', 'Resolution', 300);


%%% Plot power coeff. %%%

fig = figure('Name', 'Power coeff.', ...
    'Position', [100 100 900 350], 'Units', 'pixels');

ax1 = subplot(1,2,1);  title(ax1, title_str);
xlabel(ax1, 'Frequency (GHz)');  ylabel(ax1, 'Transmission (dB)');
hold(ax1, 'on');
ax2 = subplot(1,2,2);  title(ax2, title_str);
xlabel(ax2, 'Frequency (GHz)');  ylabel(ax2, 'Reflection (dB)');
hold(ax2, 'on');

for i = 1:numel(materialNames)
    matName = materialNames{i};
    plot(ax1, f./1e9, dBp(Tp.(matName)(ind_h,:)), 'Color', colors(i));
    plot(ax2, f./1e9, dBp(Rp.(matName)(ind_h,:)), 'Color', colors(i));
end
linkaxes([ax1, ax2], 'x');  ylim(ax1, [-40, 0]);  ylim(ax2, [-16, 0]);

legend(ax1, {"Low water content", "Medium water content", "High water content"}, ...
    'Location', 'best');
% exportgraphics(fig, '../Figures/power_RT.png', 'Resolution', 300);


%%% Compare different thickness %%%

fig = figure('Name', 'Compare thickness', ...
    'Position', [100 100 900 350], 'Units', 'pixels');

ax1 = subplot(1,2,1);  title(ax1, 'Transmission');
xlabel(ax1, 'Frequency (GHz)');  ylabel(ax1, 'Transmission (dB)');
hold(ax1, 'on');
ax2 = subplot(1,2,2);  title(ax2, 'Reflection');
xlabel(ax2, 'Frequency (GHz)');  ylabel(ax2, 'Reflection (dB)');
hold(ax2, 'on');

for i = 1:length(h)

    plot(ax1, f./1e9, dBp(Tp.lwc(i,:)), 'Color', colors(i));
    plot(ax2, f./1e9, dBp(Rp.lwc(i,:)), 'Color', colors(i));

end
linkaxes([ax1, ax2], 'x');  ylim(ax1, [-20, 0]);  ylim(ax2, [-20, 0]);

legend(ax1, {"1 cm", "2 cm", "5 cm"}, 'Location', 'best');
% exportgraphics(fig, '../Figures/power_compare.png', 'Resolution', 300);


%% Ray path

% Source & Receiver locations (m)
yt = 2e-2;  yr = -2e-2;

% Thickness of tissue (m)
h = 1e-2;

% Incidence angle (deg)
theta_i = 10;

% Refraction angles at frequencies (deg)
theta_t = atand(sind(theta_i) ./ real(sqrt(epsd.lwc - sind(theta_i)^2)));

% Ray path
y = [yt, h/2, -h/2, yr];
xinc1 = (yt - h/2) * tand(theta_i);
xinc2 = xinc1 + h .* tand(theta_t);
xr = xinc2 + (-h/2 - yr) * tand(theta_i);

% Visualize ray paths
figure('Name', 'Ray paths');
plot([0, xinc1, xinc2(1), xr(1)], y, 'k-');  hold on;
plot([0, xinc1, xinc2(end), xr(end)], y, 'k-');

% MATLAB script to calculate and plot dry peridotite solidus
% Reference: Hirschmann, M. M. (2000). Mantle solidus: Experimental constraints 
% and the effects of peridotite composition. Geochemistry, Geophysics, Geosystems, 1(10).
% 
clear; clc; close all;

% Define the pressure range (P in GPa)
P = linspace(0, 10, 100); % Pressure from 0 to 10 GPa

% Calculate temperature using the solidus equation
T = -5.104 * P.^2 + 132.899 * P + 1120.661; % Temperature in °C

% Constants for depth calculation
% Reference: Turcotte, D. L., & Schubert, G. (2014). Geodynamics (3rd ed.). Cambridge University Press.
rho = 3300; % Density of upper mantle in kg/m^3
g = 9.81;   % Gravity in m/s^2

% Depth calculation (in km)
Depth = (P * 1e9) / (rho * g) / 1000;

% Create high-quality figure
figure('Position', [100, 100, 900, 700], 'Color', 'white');

% Primary plot: Temperature vs. Pressure
yyaxis left;
plot(T, P, 'LineWidth', 3, 'Color', 'b');
ylabel('Pressure (GPa)', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'YDir', 'reverse');
ylim([0, 10]);

% Secondary plot: Depth axis
yyaxis right;
plot(T, Depth, 'LineWidth', 0.1, 'Color', 'r', 'Visible', 'off'); % Invisible plot for axis
ylabel('Depth (km)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'r');
set(gca, 'YDir', 'reverse');
ylim([0, max(Depth)]);

% Common plot settings
xlabel('Temperature (°C)', 'FontSize', 14, 'FontWeight', 'bold');
title('Dry Peridotite Solidus', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'GridAlpha', 0.3);

% Add annotation with the equation
annotation('textbox', [0.15, 0.15, 0.3, 0.1], 'String', ...
    'T = -5.104P² + 132.899P + 1120.661 (Hirschmann, M. M. , 2000)', ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
    'FontSize', 10);

% Export with high resolution
print('dry_peridotite_solidus', '-dpng', '-r300');
savefig('dry_peridotite_solidus.fig');

% Display summary statistics
fprintf('Dry Peridotite Solidus Summary:\n');
fprintf('Temperature range: %.1f°C to %.1f°C\n', min(T), max(T));
fprintf('Depth range: 0 km to %.1f km\n', max(Depth));
fprintf('Pressure range: 0 GPa to 10 GPa\n');
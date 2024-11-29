
% Constants
T = 1573;              % Temperature in Kelvin
P = 2e9;               % Pressure in Pa
Xfe = 1.0;             % Fe content (dimensionless)
Cw_values = logspace(-6, 2, 100); % Water content range in wt% (log scale from 0.0001 to 1)

% Initialize arrays to store conductivity values
sigma_Yoshino = zeros(size(Cw_values));
sigma_Karato = zeros(size(Cw_values));

% Calculate conductivity for each water content value
for i = 1:length(Cw_values)
    Cw = Cw_values(i);
    % Use Yoshino olivine conductivity subroutine
    sigma_Yoshino(i) = Yoshino_olivine_conductivity(T, Cw, P, Xfe);
    % Use Karato olivine conductivity subroutine
    sigma_Karato(i) = Karato_olivine_conductivity(T, Cw, P);
end

% Plot the results with log-log scale
figure('Units', 'inches', 'Position', [0 0 6 6.5 ], 'PaperPositionMode', 'auto');

loglog(Cw_values, sigma_Yoshino, '--r', 'LineWidth', 4); % Yoshino curve
hold on;
loglog(Cw_values, sigma_Karato, '--k', 'LineWidth', 4); % Karato curve
hold off;

% Add labels, legend, and title
xlabel('Water Content (wt%)');
ylabel('Conductivity (S/m)');
title('Olivine Conductivity via Water Content');
legend('Yoshino Olivine Conductivity', 'Wang and Karato Olivine Conductivity', 'Location', 'Best');
grid on;

% Set axis limits and ticks for better visualization
xlim([1e-4, 1e0]); % Water content from 0.01 wt% to 100 wt%
ylim([1e-3, 1e2]); % Conductivity range (adjust based on expected results)
set(gca, 'FontSize', 14); % Increase font size for better readability

% Save high-quality plot
set(gcf, 'PaperPositionMode', 'auto'); % Use screen size for printing
print('Yosino_Olivine_Conductivity', '-dpng', '-r300'); % Save as PNG with 300 dpi


% Clear workspace and initialize
clear;
clc;
% Constants for both plots
P = 15e9;                   % Pressure in Pa (15 GPa)
Xfe = 0.1;                  % Fe content for Yoshino model (dimensionless)
T_values = linspace(500, 2000, 100); % Temperature range in K
Cw_values = linspace(1e-3, 5, 100);   % Water content range in wt% for second plot

% Initialize arrays for conductivity values
sigma_Yoshino_0wt = zeros(size(T_values));
sigma_Yoshino_0_1wt = zeros(size(T_values));
sigma_Karato_0wt = zeros(size(T_values));
sigma_Karato_0_1wt = zeros(size(T_values));

% Plot 1: Conductivity vs. 1000/T
for i = 1:length(T_values)
    T = T_values(i); % Use T directly in Kelvin
    sigma_Yoshino_0wt(i) = Yoshino_wadsleyite_conductivity(T, 0, P, Xfe);    % Cw = 0 wt%
    sigma_Yoshino_0_1wt(i) = Yoshino_wadsleyite_conductivity(T, 0.1, P, Xfe); % Cw = 0.1 wt%
    sigma_Karato_0wt(i) = Karato_wadsleyite_conductivity(T, 0, P);          % Cw = 0 wt%
    sigma_Karato_0_1wt(i) = Karato_wadsleyite_conductivity(T, 0.1, P);      % Cw = 0.1 wt%
end

figure('Units', 'inches', 'Position', [0 0 6 6.5 ], 'PaperPositionMode', 'auto');
% Define RGB color for #EDB120
color_EDB120 = [237, 177, 32] / 255; % Convert from 0-255 range to 0-1

% Plot with the custom color
semilogy(1000 ./ T_values, sigma_Yoshino_0wt, '-r', 'LineWidth', 2); hold on;
semilogy(1000 ./ T_values, sigma_Yoshino_0_1wt, '--r', 'LineWidth', 2);
semilogy(1000 ./ T_values, sigma_Karato_0wt, '-', 'Color', color_EDB120, 'LineWidth', 2);
semilogy(1000 ./ T_values, sigma_Karato_0_1wt, '--', 'Color', color_EDB120, 'LineWidth', 2);
hold off;

% Customize the first plot
xlabel('1000 / T (K^{-1})');
ylabel('Conductivity (\sigma) [S/m]');
title('Wadsleyite Conductivity vs. 1000/T');
legend('Yoshino (Cw = 0 wt%)', 'Yoshino (Cw = 0.1 wt%)', ...
       'Karato (Cw = 0 wt%)', 'Karato (Cw = 0.1 wt%)', 'Location', 'Best');
grid on;
xlim([0.5, 1.0]); % Adjusted for 1000/T range
ylim([1e-4, 1e0]); % Conductivity range
set(gca, 'FontSize', 14); % Increase font size for readability

% Save first plot as high-quality image
set(gcf, 'PaperPositionMode', 'auto');
print('Yosino_wadsleyite_conductivity_vs_1000_T_4Curves', '-dpng', '-r300');

% Plot 2: Conductivity vs. Water Content
T_fixed = 1700; % Fixed temperature for second plot
sigma_Yoshino_Cw = zeros(size(Cw_values));
sigma_Karato_Cw = zeros(size(Cw_values));

for i = 1:length(Cw_values)
    sigma_Yoshino_Cw(i) = Yoshino_wadsleyite_conductivity(T_fixed, Cw_values(i), P, Xfe);
    sigma_Karato_Cw(i) = Karato_wadsleyite_conductivity(T_fixed, Cw_values(i), P);
end

% Create second plot: Conductivity vs. Water Content
figure('Units', 'inches', 'Position', [0 0 6 6.5 ], 'PaperPositionMode', 'auto');
loglog(Cw_values, sigma_Yoshino_Cw, '-r', 'LineWidth', 2); hold on;
loglog(Cw_values, sigma_Karato_Cw, '-', 'Color', color_EDB120, 'LineWidth', 2);
hold off;

% Customize the second plot
xlabel('Water Content (wt%)');
ylabel('Conductivity (\sigma) [S/m]');
title('Wadsleyite at 1700 K');
legend('Yoshino', 'Karato', 'Location', 'Best');
grid on;
xlim([1e-3, 2]);
ylim([1e-4, 10]);
set(gca, 'FontSize', 14); % Increase font size for readability

% Save second plot as high-quality image
set(gcf, 'PaperPositionMode', 'auto');
print('Yosino_wadsleyite_conductivity_vs_WaterContent', '-dpng', '-r300');


% Clear workspace and initialize
clear;
clc;
% Constants for both plots
P = 15e9;                   % Pressure in Pa (15 GPa)
Xfe = 1.0;                  % Fe content for Yoshino model (dimensionless)
T_values = linspace(500, 2000, 100); % Temperature range in K
Cw_values = linspace(1e-3, 5, 100);   % Water content range in wt% for second plot

% Initialize arrays for conductivity values
sigma_Yoshino_0wt = zeros(size(T_values));
sigma_Yoshino_0_1wt = zeros(size(T_values));
sigma_Karato_0wt = zeros(size(T_values));
sigma_Karato_0_1wt = zeros(size(T_values));

% Plot 1: Conductivity vs. 1000/T
for i = 1:length(T_values)
    T = T_values(i); % Use T directly in Kelvin
    sigma_Yoshino_0wt(i) = Yoshino_ringwoodite_conductivity(T, 0, P, Xfe);    % Cw = 0 wt%
    sigma_Yoshino_0_1wt(i) = Yoshino_ringwoodite_conductivity(T, 0.1, P, Xfe); % Cw = 0.1 wt%   
    sigma_Karato_0wt(i) = Karato_ringwoodite_conductivity(T, 0.001, P);          % Cw = 0 wt%
    sigma_Karato_0_1wt(i) = Karato_ringwoodite_conductivity(T, 0.1, P);      % Cw = 0.1 wt%
end

figure('Units', 'inches', 'Position', [0 0 6 6.5 ], 'PaperPositionMode', 'auto');
% Define RGB color for #EDB120
color_EDB120 = [237, 177, 32] / 255; % Convert from 0-255 range to 0-1

% Plot with the custom color
semilogy(1000 ./ T_values, sigma_Yoshino_0wt, '-r', 'LineWidth', 2); hold on;
semilogy(1000 ./ T_values, sigma_Yoshino_0_1wt, '--r', 'LineWidth', 2);
semilogy(1000 ./ T_values, sigma_Karato_0wt, '-', 'Color', color_EDB120, 'LineWidth', 2);
semilogy(1000 ./ T_values, sigma_Karato_0_1wt, '--', 'Color', color_EDB120, 'LineWidth', 2);
hold off;

% Customize the first plot
xlabel('1000 / T (K^{-1})');
ylabel('Conductivity (\sigma) [S/m]');
title('ringwoodite Conductivity vs. 1000/T');
legend('Yoshino (Cw = 0 wt%)', 'Yoshino (Cw = 0.1 wt%)', ...
       'Karato (Cw = 0 wt%)', 'Karato (Cw = 0.1 wt%)', 'Location', 'Best');
grid on;
xlim([0.5, 1.0]); % Adjusted for 1000/T range
ylim([1e-4, 1e0]); % Conductivity range
set(gca, 'FontSize', 14); % Increase font size for readability

% Save first plot as high-quality image
set(gcf, 'PaperPositionMode', 'auto');
print('Yosino_ringwoodite_conductivity_vs_1000_T_4Curves', '-dpng', '-r300');

% Plot 2: Conductivity vs. Water Content
T_fixed = 1700; % Fixed temperature for second plot
sigma_Yoshino_Cw = zeros(size(Cw_values));
sigma_Karato_Cw = zeros(size(Cw_values));

for i = 1:length(Cw_values)
    sigma_Yoshino_Cw(i) = Yoshino_ringwoodite_conductivity(T_fixed, Cw_values(i), P, Xfe);
    sigma_Karato_Cw(i) = Karato_ringwoodite_conductivity(T_fixed, Cw_values(i), P);
end

% Create second plot: Conductivity vs. Water Content
figure('Units', 'inches', 'Position', [0 0 6 6.5 ], 'PaperPositionMode', 'auto');
loglog(Cw_values, sigma_Yoshino_Cw, '-r', 'LineWidth', 2); hold on;
loglog(Cw_values, sigma_Karato_Cw, '-', 'Color', color_EDB120, 'LineWidth', 2);
hold off;

% Customize the second plot
xlabel('Water Content (wt%)');
ylabel('Conductivity (\sigma) [S/m]');
title('ringwoodite at 1700 K');
legend('Yoshino', 'Karato', 'Location', 'Best');
grid on;
xlim([1e-3, 2]);
ylim([1e-4, 10]);
set(gca, 'FontSize', 14); % Increase font size for readability

% Save second plot as high-quality image
set(gcf, 'PaperPositionMode', 'auto');
print('Yosino_ringwoodite_Conductivity_vs_WaterContent', '-dpng', '-r300');


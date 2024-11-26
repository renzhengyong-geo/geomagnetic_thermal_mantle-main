% Script to test the results shown in the citation:
% [1] Karato, S. I. (2011). Water distribution across the mantle transition zone
% and its implications for global material circulation. Earth and Planetary
% Science Letters, 301(3–4), 413–423. https://doi.org/10.1016/j.epsl.2010.11.038

clear;

% -------------------------------------------------------------------------
% Test the data presented in Figure 2 of the paper.
% -------------------------------------------------------------------------

% Parameters
T = 1500;                           % Temperature in Kelvin
Cw = logspace(-5, log10(1), 50);    % Water content from 10^-5 to 2 wt%
xn = length(Cw);                    % Number of water content points

% Initialize arrays for conductivities of different minerals
olivine = zeros(xn, 1);
garnet = zeros(xn, 1);
wadsleyite = zeros(xn, 1);
ringwoodite = zeros(xn, 1);
orthopyroxene = zeros(xn, 1);
clinopyroxene = zeros(xn, 1);

% -------------------------------------------------------------------------
% Compute conductivity for each mineral across water content values
% -------------------------------------------------------------------------
for i = 1:xn
    olivine(i) = Karato_olivine_conductivity(T, Cw(i), 4 * 1e9);          % 4 GPa
    garnet(i) = Karato_garnet_conductivity(T, Cw(i), 4 * 1e9);            % 4 GPa
    % orthopyroxene(i) = Karato_orthopyroxene_conductivity(T, Cw(i), 4 * 1e9); % 4 GPa
    % clinopyroxene(i) = Karato_clinopyroxene_conductivity(T, Cw(i), 4 * 1e9); % 4 GPa
    wadsleyite(i) = Karato_wadsleyite_conductivity(T, Cw(i), 15 * 1e9);   % 15 GPa
    % ringwoodite(i) = Karato_ringwoodite_conductivity(T, Cw(i), 15 * 1e9); % 15 GPa
end

% -------------------------------------------------------------------------
% High-quality plot settings
% -------------------------------------------------------------------------
figure('Units', 'inches', 'Position', [0 0 7 6], 'PaperPositionMode', 'auto');

% Plot each mineral's conductivity
plot(Cw, olivine, '-', 'LineWidth', 2); hold on;
plot(Cw, garnet, '-', 'LineWidth', 2);
% plot(Cw, orthopyroxene, '-', 'LineWidth', 2);
% plot(Cw, clinopyroxene, '-', 'LineWidth', 2);
plot(Cw, wadsleyite, '-', 'LineWidth', 2);
% plot(Cw, ringwoodite, '-', 'LineWidth', 2);
hold off;

% Set logarithmic scale for axes
set(gca, 'YScale', 'log', 'XScale', 'log', 'FontSize', 12, 'FontName', 'Helvetica');

% Axis labels and title
xlabel('C_w (wt%)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Conductivity (S/m)', 'FontSize', 14, 'FontWeight', 'bold');
title('T = 1500 K', 'FontSize', 16, 'FontWeight', 'bold');

% Define plot limits
xlim([1e-5, 1]);
ylim([1e-5, 10]);

% Add legend
% legend('Olivine', 'Garnet', 'Orthopyroxene', 'Clinopyroxene', 'Wadsleyite', 'Ringwoodite', ...
%        'Location', 'best', 'FontSize', 12);
legend('Olivine', 'Garnet', 'Wadsleyite', 'Location', 'best', 'FontSize', 12);
% Improve plot appearance
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'TickLength', [0.02 0.02]);
box on;

% -------------------------------------------------------------------------
% Export the plot as a high-resolution image
% -------------------------------------------------------------------------
print(gcf, 'Conductivity_vs_WaterContent', '-dpng', '-r300'); % 300 DPI for publication

% -------------------------------------------------------------------------
% Test Olivine data in Figure 2b of Wang et al.(2006)
% Wang, D., Mookherjee, M., Xu, Y., & Karato, S. I. (2006). The effect of 
% water on the electrical conductivity of olivine. 
% Nature, 443(7114), 977–980. https://doi.org/10.1038/nature05256
% -------------------------------------------------------------------------
% Clear workspace and initialize
clear;
clc;

% Define constants
P = 3e9; % Pressure in Pascal (constant for upper mantle)
T_values = [1273, 1173, 1073, 973, 873]; % Temperatures in Kelvin
Cw = logspace(log10(1e-3), log10(1), 50); % Water content in wt%
colors = lines(length(T_values)); % Generate colors for the lines

% Load experimental data
data = readmatrix('wang_nature_2006_fig2b.csv');
exp_Cw = data(:, 1); % Experimental water content
exp_sigma = data(:, 2); % Experimental conductivity

% Create figure
figure('Color', 'white', 'Units', 'inches', 'Position', [1, 1, 7, 7]);
hold on;

% Plot theoretical curves
for i = 1:length(T_values)
    T = T_values(i); % Current temperature
    sigma = zeros(size(Cw)); % Preallocate conductivity array
    for j = 1:length(Cw)
        sigma(j) = Karato_olivine_conductivity(T, Cw(j), P); % Compute conductivity
    end
    % Plot conductivity curve
    plot(Cw, sigma, '-', 'LineWidth', 2, 'DisplayName', sprintf('T = %d K', T), 'Color', colors(i, :));
end

% Plot experimental data
scatter(exp_Cw, exp_sigma, 50, 'k', 'filled', 'DisplayName', 'Exp. Data (Wang et al., 2006)');

% Customize plot
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 12, 'LineWidth', 1.5);
xlim([1e-3, 1]);
ylim([5e-5, 2]);
xlabel('Water Content (wt%)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Conductivity (S/m)', 'FontSize', 14, 'FontWeight', 'bold');
title('Olivine Conductivity vs. Water Content', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
grid on;
box on;

% Export the plot as a high-resolution image
print(gcf, 'Karato_Wang_2006_Olivine_Conductivity_vs_WaterContent', '-dpng', '-r300');


% -------------------------------------------------------------------------
% Test Opx and Garnet (also olivine) data in Figure 3 of the following paper
% Dai, L., & Karato, S. ichiro. (2009). Electrical conductivity of orthopyroxene: 
% Implications for the water content of the asthenosphere. Proceedings of the Japan 
% Academy Series B: Physical and Biological Sciences, 85(10), 466–475. 
% https://doi.org/10.2183/pjab.85.466
% Dai, L., & Karato, S. ichiro. (2009). Electrical conductivity of pyrope-rich garnet 
% at high temperature and high pressure. Physics of the Earth and Planetary Interiors,
% 176(1–2), 83–88. https://doi.org/10.1016/j.pepi.2009.04.002
% -------------------------------------------------------------------------
% Clear workspace and initialize
clear;
clc;

% Define constants
P = 8e9; % Pressure in Pascal (constant)
Cw_values = [0.042,  0]; % Water content (wt%)
T_values = linspace(833, 1666, 100); % T range corresponding to 10000/T from 6 to 12 (K)

% Compute 10000/T
x_axis = 10000 ./ T_values;

% Load experimental data
data = readmatrix('Dai_Karato_PJAS_2009_fig3.csv');
exp_x = data(:, 1); % Experimental 10000/T
exp_sigma = data(:, 2); % Experimental conductivity

% Preallocate conductivities
sigma_orthopyroxene = zeros(length(Cw_values), length(T_values));
sigma_garnet = zeros(length(Cw_values), length(T_values));
sigma_olivine = zeros(length(Cw_values), length(T_values));

% Compute theoretical conductivities for each water content
for i = 1:length(Cw_values)
    Cw = Cw_values(i);
    for j = 1:length(T_values)
        T = T_values(j);
        sigma_orthopyroxene(i, j) = Karato_orthopyroxene_conductivity(T, Cw, P);
        sigma_garnet(i, j) = Karato_garnet_conductivity(T, Cw, P);
        sigma_olivine(i, j) = Karato_olivine_conductivity(T, Cw, P);
    end
end

% Create the plot
figure('Color', 'white', 'Units', 'inches', 'Position', [1, 1, 7, 7]);

% Plot theoretical curves for orthopyroxene
plot(x_axis, sigma_garnet(1, :), '-.', 'LineWidth', 2, 'DisplayName', 'Gt: C_w = 0.042 wt%');
hold on;
plot(x_axis, sigma_orthopyroxene(1, :), '-', 'LineWidth', 2, 'DisplayName', 'Opx: C_w = 0.042 wt%');
plot(x_axis, sigma_olivine(1, :), '-o', 'LineWidth', 2, 'DisplayName', 'Ol: C_w = 0.042 wt%');
plot(x_axis, sigma_garnet(2, :), ':', 'LineWidth', 2, 'DisplayName', 'Gt: C_w = 0 wt%');
plot(x_axis, sigma_orthopyroxene(2, :), '--', 'LineWidth', 2, 'DisplayName', 'Opx: C_w = 0 wt%');
plot(x_axis, sigma_olivine(2, :), '--o', 'LineWidth', 2, 'DisplayName', 'Ol: C_w = 0 wt%');

% Plot experimental data
scatter(exp_x, exp_sigma, 50, 'k', 'filled', 'DisplayName', 'Exp. Data (Dai & Karato, 2009)');

% Customize plot
set(gca, 'YScale', 'log', 'FontSize', 12, 'LineWidth', 1.5);
xlim([6, 12]);
ylim([1e-8, 1]);
xlabel('10000 / T (K^{-1})', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Conductivity (S/m)', 'FontSize', 14, 'FontWeight', 'bold');
title('Conductivity vs. 10000/T (Opx, Gt, Ol)', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
grid on;
box on;

% Export the plot as a high-resolution image
print(gcf, 'Dai_Karato_2009_Opx_Gt_Ol_Conductivity_vs_10000_T', '-dpng', '-r300');


% -------------------------------------------------------------------------
% Test wadsleyite data in Figure 6 of the following paper
% Dai, L., & Karato, S. ichiro. (2009). Electrical conductivity of wadsleyite 
% at high temperatures and high pressures. Earth and Planetary Science Letters, 
% 287(1–2), 277–283. https://doi.org/10.1016/j.epsl.2009.08.012
% -------------------------------------------------------------------------
% Clear workspace and initialize
clear;
clc;

% Define constants
P = 15e9; % Pressure in Pascal (for wadsleyite, typical for mantle transition zone)
Cw_values = logspace(-3, 0, 100); % Water content range from 10^-3 to 1 wt%
T_values = [1273, 1173, 1073, 973, 873]; % Temperatures in Kelvin

% Preallocate for conductivities
sigma_wadsleyite = zeros(length(T_values), length(Cw_values));

% Compute conductivity for each temperature and water content
for i = 1:length(T_values)
    T = T_values(i);
    for j = 1:length(Cw_values)
        sigma_wadsleyite(i, j) = Karato_wadsleyite_conductivity(T, Cw_values(j), P);
    end
end

% Load experimental data
data = readmatrix('Dai_Karato_EPSL_2009_fig6.csv');
exp_Cw = data(:, 1); % Experimental water content (wt%)
exp_sigma = data(:, 2); % Experimental conductivity (S/m)

% Create the plot
figure('Color', 'white', 'Units', 'inches', 'Position', [1, 1, 7, 7]);

% Plot theoretical curves
hold on;
colors = lines(length(T_values)); % Generate unique colors
for i = 1:length(T_values)
    plot(Cw_values, sigma_wadsleyite(i, :), '-', 'LineWidth', 2, 'Color', colors(i, :), ...
         'DisplayName', sprintf('T = %d K', T_values(i)));
end

% Plot experimental data
scatter(exp_Cw, exp_sigma, 50, 'k', 'filled', 'DisplayName', 'Exp. Data (Dai & Karato, EPSL, 2009)');

% Customize plot
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 12, 'LineWidth', 1.5);
xlim([1e-3, 1]);
ylim([1e-5, 1]);
xlabel('Water Content (wt%)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Conductivity (S/m)', 'FontSize', 14, 'FontWeight', 'bold');
title('Conductivity vs. Water Content (Wadsleyite)', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
grid on;
box on;

% Export the plot as a high-resolution image
print(gcf, 'Wadsleyite_Conductivity_vs_WaterContent', '-dpng', '-r300');

% -------------------------------------------------------------------------
% Test ringwoodite data in Figure 1 of the following paper
% Huang, X., Xu, Y., & Karato, S. (2005). Water content in the transition zone 
% from electrical conductivity of wadsleyite and ringwoodite. 
% Nature, 434(7034), 746–749. https://doi.org/10.1038/nature03426
% -------------------------------------------------------------------------
% Clear workspace and initialize
clear;
clc;

% Define constants and ranges
T_values = linspace(833, 1666, 100); % Temperature range in K
water_contents = logspace(-2, 1, 100); % Water content range in wt%
P = 15 * 1e9; % Pressure in Pa

% Define water contents for subplots 1 and 3
Cw_values = [0.01, 0.1, 1, 10]; 

% Create figure
figure('Color', 'white', 'Units', 'inches', 'Position', [1, 1, 12, 8]);

% Subplot 1: Wadsleyite conductivity vs 10000/T
subplot(2, 2, 1);
hold on;
for Cw = Cw_values
    conductivity = arrayfun(@(T) Karato_wadsleyite_conductivity(T, Cw, P), T_values);
    plot(10000 ./ T_values, conductivity, 'LineWidth', 2, 'DisplayName', sprintf('Cw = %.2f wt%%', Cw));
end
hold off;
set(gca, 'YScale', 'log', 'FontSize', 12);
xlabel('10000/T (K^{-1})', 'FontSize', 14);
ylabel('Conductivity (S/m)', 'FontSize', 14);
title('Wadsleyite Conductivity vs 10000/T', 'FontSize', 16);
legend('Location', 'best', 'FontSize', 10);
xlim([6, 12]);
ylim([1e-4, 1]);
grid on;

% Subplot 2: Wadsleyite conductivity vs water content
subplot(2, 2, 2);
hold on;
for T = [873, 1073, 1273, 1473]
    conductivity = arrayfun(@(Cw) Karato_wadsleyite_conductivity(T, Cw, P), water_contents);
    plot(water_contents, conductivity, 'LineWidth', 2, 'DisplayName', sprintf('T = %d K', T));
end
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 12);
xlabel('Water Content (wt%)', 'FontSize', 14);
ylabel('Conductivity (S/m)', 'FontSize', 14);
title('Wadsleyite Conductivity vs Water Content', 'FontSize', 16);
legend('Location', 'best', 'FontSize', 10);
grid on;
ylim([1e-4, 1]);

% Subplot 3: Ringwoodite conductivity vs 10000/T
subplot(2, 2, 3);
hold on;
for Cw = Cw_values
    conductivity = arrayfun(@(T) Karato_ringwoodite_conductivity(T, Cw, P), T_values);
    plot(10000 ./ T_values, conductivity, 'LineWidth', 2, 'DisplayName', sprintf('Cw = %.2f wt%%', Cw));
end
hold off;
set(gca, 'YScale', 'log', 'FontSize', 12);
xlabel('10000/T (K^{-1})', 'FontSize', 14);
ylabel('Conductivity (S/m)', 'FontSize', 14);
title('Ringwoodite Conductivity vs 10000/T', 'FontSize', 16);
legend('Location', 'best', 'FontSize', 10);
xlim([6, 12]);
ylim([1e-4, 1]);
grid on;

% Subplot 4: Ringwoodite conductivity vs water content
subplot(2, 2, 4);
hold on;
for T = [873, 1073, 1273, 1473]
    conductivity = arrayfun(@(Cw) Karato_ringwoodite_conductivity(T, Cw, P), water_contents);
    plot(water_contents, conductivity, 'LineWidth', 2, 'DisplayName', sprintf('T = %d K', T));
end
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 12);
xlabel('Water Content (wt%)', 'FontSize', 14);
ylabel('Conductivity (S/m)', 'FontSize', 14);
title('Ringwoodite Conductivity vs Water Content', 'FontSize', 16);
legend('Location', 'best', 'FontSize', 10);
grid on;
ylim([1e-4, 1]);

% Adjust layout
sgtitle('Conductivity Analysis for Wadsleyite and Ringwoodite', 'FontSize', 18);


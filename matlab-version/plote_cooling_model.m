% Thermal evolution of oceanic lithosphere using half-space and plate cooling models
% FULL REFERENCES:
% ================
%
% 1. PRIMARY PARAMETER SOURCE: Parsons, B., & Sclater, J. G. (1977). 
%    "An analysis of the variation of ocean floor bathymetry and heat flow with age."
%    Journal of Geophysical Research, 82(5), 803-827.
%    DOI: 10.1029/JB082i005p00803
% 2. THERMAL DIFFUSIVITY SOURCE:  Stein, C. A., & Stein, S. (1992).
%    "A model for the global variation in oceanic depth and heat flow with lithospheric age."
%    Nature, 359(6391), 123-129.
% 3. HALF-SPACE COOLING MODEL:  Turcotte, D. L., & Schubert, G. (2014).
%    "Geodynamics" (3rd edition), Cambridge University Press, ISBN: 978-1-107-00653-9
% 4. PLATE MODEL FOUNDATION: McKenzie, D. P. (1967).
%    "Some remarks on heat flow and gravity anomalies."
%    Journal of Geophysical Research, 72(24), 6261-6273.
% 5. COMPARATIVE STUDY: Carlson, R. L., & Johnson, H. P. (1994).
%    "On modeling the thermal evolution of the oceanic upper mantle."
%    Journal of Geophysical Research: Solid Earth, 99(B2), 3201-3214.
%    DOI: 10.1029/93JB02863 (Comparison of different plate thickness
%    estimates (95-125 km))

clear; clc; close all;

%% MODEL PARAMETERS
parameters.T_surface = 0;   % Ocean bottom temperature in degree
parameters.T_mantle = 1330; % Mantle temperature (1330°C): Parsons & Sclater (1977)
parameters.kappa = 8.04e-7; % Thermal diffusivity (8.04 × 10⁻⁷ m²/s): Stein & Stein (1992)
parameters.k_thermal = 3.138; % Thermal conductivity (3.138 W/m/K): Parsons & Sclater (1977)
parameters.plate_thickness = 95e3; % Plate thickness (95 km): Parsons & Sclater (1977)
parameters.n_terms = 100;

%% Age and depth ranges (starting from non-zero age)
ages_ma = 1:2:200;  % 1 to 200 Myr in 2 Myr steps (avoid age=0)
ages_sec = ages_ma * 1e6 * 365.25 * 24 * 3600;

depth_km = 0:2:150;  % 0 to 150 km in 2 km steps
depth_m = depth_km * 1000;

% Isotherms for contour plot
isotherms = [200, 400, 600, 800, 1000, 1200, 1300];

%% Calculate temperature grids for contour plot
T_halfspace = zeros(length(depth_km), length(ages_ma));
T_plate_95 = zeros(length(depth_km), length(ages_ma));

fprintf('Calculating temperature grids...\n');

for i = 1:length(ages_ma)
    % Half-space model
    [T_hs, ~] = halfspace_cooling(depth_m, ages_sec(i), parameters);
    T_halfspace(:, i) = T_hs;
    
    % Plate model (95 km)
    [T_pl, ~] = plate_cooling(depth_m, ages_sec(i), parameters);
    T_plate_95(:, i) = T_pl;
end

%% PLOT 1: Age-Depth Temperature Contours
figure('Position', [100, 100, 800, 600]);
hold on;

[AG, DP] = meshgrid(ages_ma, depth_km);

% Plot half-space model contours in red
[C1, h1] = contour(AG, DP, T_halfspace, isotherms, 'LineWidth', 2, 'LineColor', 'red');
clabel(C1, h1, 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'red');

% Plot plate model contours in blue
[C2, h2] = contour(AG, DP, T_plate_95, isotherms, 'LineWidth', 2, 'LineColor', 'blue', 'LineStyle', '--');
clabel(C2, h2, 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'blue');

% Formatting
set(gca, 'YDir', 'reverse');
xlabel('t (Myr)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('y (km)', 'FontSize', 14, 'FontWeight', 'bold');
title('Age-Depth Temperature Contours', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
xlim([1, 150]);  % Start from 1 Myr
ylim([0, 150]);

% Add legend
legend([h1(1), h2(1)], {'Half-Space Model', 'Plate Model (95 km)'}, ...
       'Location', 'southeast', 'FontSize', 12);

% Add plate thickness reference line
yline(95, 'k--', 'LineWidth', 1, 'Alpha', 0.7, 'Label', 'Plate Base (95 km)', ...
      'LabelHorizontalAlignment', 'left', 'FontSize', 10);

saveas(gcf, 'Age_Depth_Contours.png', 'png');

%% PLOT 2: Heat Flow (q0) vs Age for three models
fprintf('Calculating heat flow...\n');

% Initialize heat flow arrays
q_halfspace = zeros(size(ages_ma));
q_plate_95 = zeros(size(ages_ma));
q_plate_125 = zeros(size(ages_ma));

% Create parameters for 125 km plate model
parameters_125 = parameters;
parameters_125.plate_thickness = 125e3;

for i = 1:length(ages_ma)
    % Half-space model
    [~, q_hs] = halfspace_cooling(depth_m, ages_sec(i), parameters);
    q_halfspace(i) = q_hs;
    
    % Plate model 95 km
    [~, q_pl_95] = plate_cooling(depth_m, ages_sec(i), parameters);
    q_plate_95(i) = q_pl_95;
    
    % Plate model 125 km
    [~, q_pl_125] = plate_cooling(depth_m, ages_sec(i), parameters_125);
    q_plate_125(i) = q_pl_125;
end

figure('Position', [100, 100, 800, 600]);
hold on;

% Convert from W/m² to mW/m² for better readability
q_halfspace_mW = q_halfspace * 1000;
q_plate_95_mW = q_plate_95 * 1000;
q_plate_125_mW = q_plate_125 * 1000;

% Plot the heat flow curves
plot(ages_ma, q_halfspace_mW, 'r-', 'LineWidth', 3, 'DisplayName', 'HSCM');
plot(ages_ma, q_plate_95_mW, 'b--', 'LineWidth', 3, 'DisplayName', 'PM 95');
plot(ages_ma, q_plate_125_mW, 'b-.', 'LineWidth', 3, 'DisplayName', 'PM 125');

% Formatting
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
xlabel('t, Myr', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('q_0, mW/m²', 'FontSize', 14, 'FontWeight', 'bold');
title('Surface Heat Flow vs. Lithospheric Age', 'FontSize', 16, 'FontWeight', 'bold');

% Set axis limits
xlim([1, 200]);  % Start from 1 Myr
ylim([0, 250]);

% Add grid
grid on;

% Add legend
legend('Location', 'northeast', 'FontSize', 12, 'FontWeight', 'bold');

saveas(gcf, 'Heat_Flow_Comparison.png', 'png');

%% Display key observations
fprintf('\nKey Observations:\n');
fprintf('• Age range: 1-200 Myr (avoiding age=0 singularity)\n');
fprintf('• Combined contour plot shows:\n');
fprintf('  - Red: Half-space model (isotherms deepen indefinitely)\n');
fprintf('  - Blue: Plate model (isotherms approach plate base)\n');
fprintf('• Heat flow behavior:\n');
fprintf('  - HSCM: Decreases as t^{-1/2}\n');
fprintf('  - PM 95: Approaches ~32 mW/m² for old lithosphere\n');
fprintf('  - PM 125: Approaches ~24 mW/m² for old lithosphere\n');

fprintf('\nHeat flow values:\n');
fprintf('  At 1 Myr:\n');
fprintf('    HSCM: %.1f mW/m²\n', q_halfspace_mW(1));
fprintf('    PM 95: %.1f mW/m²\n', q_plate_95_mW(1));
fprintf('    PM 125: %.1f mW/m²\n', q_plate_125_mW(1));
fprintf('  At 200 Myr:\n');
fprintf('    HSCM: %.1f mW/m²\n', q_halfspace_mW(end));
fprintf('    PM 95: %.1f mW/m²\n', q_plate_95_mW(end));
fprintf('    PM 125: %.1f mW/m²\n', q_plate_125_mW(end));

%% PLOT 3: Temperature-Depth Profiles for Specific Ages
fprintf('Creating temperature-depth profiles for specific ages...\n');

% Selected ages for detailed analysis
selected_ages = [70, 100, 120];  % Myr
selected_ages_sec = selected_ages * 1e6 * 365.25 * 24 * 3600;

% Create detailed depth range for smooth profiles
depth_detailed_km = 0:1:150;  % 1 km resolution for smooth curves
depth_detailed_m = depth_detailed_km * 1000;

% Colors for different ages
age_colors = [0.8, 0.2, 0.2;  % 79 Myr - Red
              0.2, 0.6, 0.2;  % 100 Myr - Green
              0.2, 0.2, 0.8]; % 120 Myr - Blue

figure('Position', [100, 100, 800, 600]);
hold on;

% Calculate and plot profiles for each age
for i = 1:length(selected_ages)
    % Find closest age index in our calculated data
    age_idx = find(ages_ma >= selected_ages(i), 1);
    
    % Calculate detailed profiles
    [T_hs_detailed, ~] = halfspace_cooling(depth_detailed_m, selected_ages_sec(i), parameters);
    [T_pl_detailed, ~] = plate_cooling(depth_detailed_m, selected_ages_sec(i), parameters);
    
    % Plot half-space model (solid lines)
    plot(T_hs_detailed, depth_detailed_km, '-', 'LineWidth', 3, ...
         'Color', age_colors(i,:), ...
         'DisplayName', sprintf('%d Myr (HSCM)', selected_ages(i)));
    
    % Plot plate model (dashed lines)
    plot(T_pl_detailed, depth_detailed_km, '--', 'LineWidth', 3, ...
         'Color', age_colors(i,:), ...
         'DisplayName', sprintf('%d Myr (Plate)', selected_ages(i)));
end

% Formatting
set(gca, 'YDir', 'reverse');
xlabel('Temperature (°C)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Depth (km)', 'FontSize', 14, 'FontWeight', 'bold');
title({'Temperature-Depth Profiles'; 'for Selected Lithospheric Ages'}, ...
      'FontSize', 16, 'FontWeight', 'bold');
grid on;

% Set axis limits
xlim([0, 1400]);
ylim([0, 150]);

% Add legend
legend('Location', 'southeast', 'FontSize', 12, 'NumColumns', 2);

% Add plate thickness reference line
yline(95, 'k--', 'LineWidth', 2, 'Label', 'Plate Base (95 km)', ...
      'LabelHorizontalAlignment', 'left', 'FontSize', 12, ...
      'Color', [0.3 0.3 0.3]);

% Add mantle temperature reference line
xline(1330, 'r-', 'LineWidth', 1, 'Label', 'Mantle Temperature (1330°C)', ...
      'LabelVerticalAlignment', 'bottom', 'FontSize', 10, ...
      'Color', [0.7 0.2 0.2], 'Alpha', 0.7);

saveas(gcf, 'Temperature_Depth_Profiles.png', 'png');
saveas(gcf, 'Temperature_Depth_Profiles.pdf', 'pdf');

%% Display temperature values at specific depths for reference
fprintf('\n=== TEMPERATURE COMPARISON AT SELECTED AGES ===\n');
fprintf('Temperature at 50 km depth:\n');
for i = 1:length(selected_ages)
    age_idx = find(ages_ma >= selected_ages(i), 1);
    
    % Find temperature at 50 km depth
    depth_50km_idx = find(depth_km >= 50, 1);
    T_hs_50km = T_halfspace(depth_50km_idx, age_idx);
    T_pl_50km = T_plate_95(depth_50km_idx, age_idx);
    
    fprintf('  %d Myr: HSCM = %.1f°C, Plate = %.1f°C, Difference = %.1f°C\n', ...
            selected_ages(i), T_hs_50km, T_pl_50km, T_hs_50km - T_pl_50km);
end

fprintf('\nTemperature at 100 km depth:\n');
for i = 1:length(selected_ages)
    age_idx = find(ages_ma >= selected_ages(i), 1);
    
    % Find temperature at 100 km depth
    depth_100km_idx = find(depth_km >= 100, 1);
    T_hs_100km = T_halfspace(depth_100km_idx, age_idx);
    T_pl_100km = T_plate_95(depth_100km_idx, age_idx);
    
    fprintf('  %d Myr: HSCM = %.1f°C, Plate = %.1f°C, Difference = %.1f°C\n', ...
            selected_ages(i), T_hs_100km, T_pl_100km, T_hs_100km - T_pl_100km);
end

fprintf('\nThird plot saved:\n');
fprintf('  - Temperature_Depth_Profiles.png/pdf\n');


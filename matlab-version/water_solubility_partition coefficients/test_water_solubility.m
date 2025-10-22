% TEST 1,
% Example from abstract and Table 1
% Kohlstedt, Kepper etc., CMP, 1996
P = 13;   % Pressure in GPa (~60 km depth)
T = 1100; % Temperature in °C (typical mantle geotherm)
[C_OH, C_H2O_ppm]=olivine_water_solubility(P, T);
% Display results
fprintf('Water Solubility in Enstatite - Kohlstedt et al. (1996) Model\n');
fprintf('Conditions: P = %.1f GPa, T = %.0f°C\n', P, T);
fprintf('----------------------------------------\n');
fprintf('olivine solubility: %.0f H/106Si\n', C_OH);
fprintf('olivine solubility: %.0f ppm H2O\n', C_H2O_ppm);
fprintf('----------------------------------------\n');


% TEST 2
% Example: Calculate water solubility in enstatite for typical mantle conditions
% Corresponding to approximately 60 km depth
P = 4;   % Pressure in GPa (~60 km depth)
T = 1000; % Temperature in °C (typical mantle geotherm)

% Calculate water solubility using Mierdel et al. (2007) model
[c_total, c_Al_free, c_Al] = enstatite_water_solubility(P, T);

% Display results
fprintf('Water Solubility in Enstatite - Mierdel et al. (2007) Model\n');
fprintf('Conditions: P = %.1f GPa, T = %.0f°C\n', P, T);
fprintf('----------------------------------------\n');
fprintf('Al-free enstatite solubility:    %.0f ppm H2O\n', c_Al_free);
fprintf('Al-coupled additional solubility: %.0f ppm H2O\n', c_Al);
fprintf('Total aluminous enstatite:       %.0f ppm H2O\n', c_total);
fprintf('----------------------------------------\n');
fprintf('Note: Al-coupled mechanism dominates at low P-T conditions\n');
fprintf('      Al-free mechanism becomes more important at high P-T\n');

% TEST 3
% Testing code for cpx_water_solubility function - OPTION A
% Generate water solubility vs pressure plot for chromian diopside

fprintf('=== Water Solubility in Chromian Diopside ===\n');
fprintf('Model: General solubility law (Bromiley et al. 2004) - OPTION A\n');
fprintf('Parameters for natural Cr-diopside sample DI-1\n\n');

% Create pressure range from 1 to 7 GPa (typical upper mantle range)
P_range = 1:0.1:7; % GPa
T_fixed = 1000;    % °C - fixed temperature

% Calculate water solubility for each pressure
water_content = zeros(size(P_range));
for i = 1:length(P_range)
    water_content(i) = cpx_water_solubility(P_range(i), T_fixed);
end

% Create the figure
figure('Position', [100, 100, 800, 600]);
plot(P_range, water_content, 'r-', 'LineWidth', 2);
grid on;

% Add labels and title
xlabel('Pressure (GPa)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Water Content (ppm H₂O)', 'FontSize', 12, 'FontWeight', 'bold');
title('Water Solubility in Chromian Diopside vs Pressure', 'FontSize', 14, 'FontWeight', 'bold');

% Add annotation with model parameters
annotation('textbox', [0.15, 0.7, 0.3, 0.2], 'String', ...
    sprintf('Model Parameters (Option A):\nT = %d°C\nA = %.2f ppm/bar^{0.5}\nn = %.1f\nΔV = %.2f cm³/mol\n\nIncorporation:\nCr³⁺ + H⁺ ↔ Si⁴⁺', ...
    T_fixed, 2.15, 0.5, 7.43), ...
    'BackgroundColor', 'white', 'FontSize', 10, 'EdgeColor', 'black');

% Add citation
annotation('textbox', [0.15, 0.15, 0.7, 0.1], 'String', ...
    'Equation and Parameters: Bromiley et al. (2004) American Mineralogist, 89, 941-949', ...
    'BackgroundColor', 'white', 'FontSize', 9, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');

% Set axis properties
xlim([1, 7]);
set(gca, 'FontSize', 11);

% Display some sample values
fprintf('Sample Calculations at T = %d°C:\n', T_fixed);
fprintf('P (GPa)   Water Content (ppm)\n');
fprintf('-----------------------------\n');
sample_indices = [1, 11, 21, 31, 41, 51, 61]; % Every 1 GPa
for i = sample_indices
    if i <= length(P_range)
        fprintf('  %.1f         %.0f\n', P_range(i), water_content(i));
    end
end

% Calculate and display the range
fprintf('\nWater solubility range: %.0f - %.0f ppm H₂O\n', ...
    min(water_content), max(water_content));

% Scientific interpretation
fprintf('\nScientific Interpretation:\n');
fprintf('• n = 0.5 indicates isolated OH groups\n');
fprintf('• Dominant mechanism: Cr³⁺ + H⁺ ↔ Si⁴⁺\n');
fprintf('• Model valid for natural chromian diopside compositions\n');
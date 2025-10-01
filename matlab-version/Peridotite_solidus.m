% Step 1: Load the CSV data (skip the header row)
data = readtable('Sarafian_Science_2017.csv', 'HeaderLines', 1);  % Skip the first row as it's the header

% Step 2: Extract the relevant columns
pressure = data{:, 1};  % First column: Pressure in GPa
dry_solidus_sarafian = data{:, 2};  % Second column: Dry solidus from Sarafian 2017
solidus_200ppm_water = data{:, 3};  % Third column: Solidus of peridotite with 200 ppm water
solidus_450ppm_water = data{:, 4};  % Fourth column: Solidus of peridotite with 450 ppm water
dry_solidus_hirschmann_2000 = data{:, 5};  % Fifth column: Dry solidus from Hirschmann et al. (2000)
solidus_200ppm_water_hirschmann_2009 = data{:, 6};  % Sixth column: Solidus with 200 ppm water from Hirschmann et al. (2009)
adiabat_1350 = data{:, 7};  % Seventh column: Adiabat mantle profile with a potential temperature of 1350°C
adiabat_1410 = data{:, 8};  % Eighth column: Adiabat mantle profile with a potential temperature of 1410°C

% Step 3: Create the plot (with a larger figure size)
figure('Position', [100, 100, 800, 700]);  % Increase figure size to 12x8 inches
hold on;

% Plotting the curves with Nature-friendly colors
plot(pressure, dry_solidus_sarafian, 'r-', 'LineWidth', 2, 'DisplayName', 'Dry Solidus (Sarafian 2017)');
plot(pressure, solidus_200ppm_water, 'r--', 'LineWidth', 2, 'DisplayName', 'Solidus (200 ppm H2O, Sarafian 2017)');  % Blue
plot(pressure, solidus_450ppm_water, 'r-.', 'LineWidth', 2, 'DisplayName', 'Solidus (450 ppm H2O, Sarafian 2017)');  % Green
plot(pressure, dry_solidus_hirschmann_2000, 'k-', 'LineWidth', 2, 'DisplayName', 'Dry Solidus (Hirschmann 2000)');  % Orange
plot(pressure, solidus_200ppm_water_hirschmann_2009, 'k--', 'LineWidth', 2, 'DisplayName', 'Solidus (200 ppm H2O, Hirschmann 2009)');  % Purple
plot(pressure, adiabat_1350, 'k-', "Marker","<", 'LineWidth', 2, 'DisplayName', 'Adiabat (1350°C)');  % Teal
plot(pressure, adiabat_1410, 'k--', "Marker","<", 'LineWidth', 2, 'DisplayName', 'Adiabat (1410°C)');  % Yellow

% Step 4: Customize the plot
xlabel('Pressure (GPa)', 'FontSize', 14);
ylabel('Temperature (°C)', 'FontSize', 14);
title('Pressure vs Temperature Solidus and Adiabat Curves', 'FontSize', 16);

% Set axis limits
xlim([0 3]);  % Set x-axis from 0 to 3 GPa
ylim([900 1600]);  % Set y-axis from 900°C to 1600°C

% Add a legend and place it outside the plot to avoid overlap
legend('show', 'Location', 'best', 'FontSize', 12);

% Set the figure properties for high-quality publication
set(gca, 'FontSize', 12);
set(gca, 'LineWidth', 1.5);
box on;  % Show the box around the plot
grid on;  % Display the grid

% Step 5: Save the figure in high-quality format for publication
print('Sarafian_Science_2017_Plot', '-dpng', '-r300');  % Save as PNG with 300 DPI (high resolution)

hold off;

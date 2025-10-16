% MATLAB code to read the CSV file and generate the plot, skipping the first row
% Make sure your CSV file is in the correct directory or provide the full path.

% Read the data from the CSV file, skipping the first row
data = readtable('Dasgupta_Hirschmann_Nature_2006_Carbonatue.csv', 'HeaderLines', 1);

% Extract the columns (temperature and pressure)
temperature = data{:, 1}; % First column is temperature in °C
pressure = data{:, 2};     % Second column is pressure in GPa

% Convert pressure to depth using the relationship: depth (km) ≈ pressure (GPa) × 30
% This is a common approximation for the upper mantle (density ~3300 kg/m³)
depth_km = pressure * 30;

% Create the plot with two y-axes
figure;

% Primary y-axis for pressure
yyaxis left
plot(temperature, pressure, 'LineWidth', 2, 'Color', 'b'); % Line plot with thicker blue line
ylabel('Pressure (GPa)', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'YDir', 'reverse');
ylim([1, 11]);      % Set y-axis from 1 GPa to 11 GPa

% Secondary y-axis for depth
yyaxis right
plot(temperature, depth_km, 'LineWidth', 2, 'Color', 'r', 'Visible', 'off'); % Invisible line for depth
ylabel('Depth (km)', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'YDir', 'reverse'); % Depth also increases downward
ylim([1*30, 11*30]); % Set depth axis from 30 km to 330 km

% Common x-axis settings
xlabel('Temperature (°C)', 'FontSize', 14, 'FontWeight', 'bold');
xlim([900, 1600]);  % Set x-axis from 900°C to 1600°C

% Title for the plot
title('25000ppm Carbonature peridotite solidus', 'FontSize', 16, 'FontWeight', 'bold');

% Display grid for better readability
grid on;

% Increase the font size for the ticks
set(gca, 'FontSize', 12);

% Save the figure with the new name
saveas(gcf, 'Carbonature_25000_peridotite_solidus_with_depth.png', 'png');

% Optional: To save in PDF format with high quality
% saveas(gcf, 'Carbonature_25000_peridotite_solidus_with_depth.pdf', 'pdf');
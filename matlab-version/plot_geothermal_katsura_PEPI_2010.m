% Define the filename
filename = 'geothermal_katsura_PEPI_2010.csv'; % Replace with your file name

% Read the CSV file, skipping the first row
opts = detectImportOptions(filename);  % Automatically detect the file structure
opts.DataLines = [2, Inf];  % Start reading from the second row
data = readtable(filename, opts);  % Read the table

% Extract depth and temperature columns
depth = data.Depth;  % Adjust if your column names are different
temperature = data.T;  % Adjust if your column names are different

% Create the T-Depth plot
figure;

% Plot with customizations for high quality
plot(depth, temperature, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]);  % Line color can be adjusted
xlabel('Depth (km)', 'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');
ylabel('Temperature (K)', 'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');
xlim([0 3000]);  % Range for depth
ylim([1000 3500]);  % Range for temperature

% Set axis properties
set(gca, 'XScale', 'linear', 'YScale', 'linear', 'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'normal');
title('Temperature vs Depth', 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');
grid on;

% Increase figure size and set resolution for high-quality output
set(gcf, 'Position', [100, 100, 800, 600]);  % Set figure size (width x height)
set(gca, 'LineWidth', 1.5);  % Increase axis line width

% Export the plot to a high-quality PNG file
exportgraphics(gcf, 'geothermal_katsura_PEPI_2010.png', 'Resolution', 300);


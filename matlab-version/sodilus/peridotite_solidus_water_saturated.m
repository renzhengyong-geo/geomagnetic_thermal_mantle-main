% MATLAB script to read and plot Kawamoto water-saturated peridotite solidus data
% Reference: Kawamoto, J. R. Holloway, Melting temperature and partial melt chemistry
% to H2O-­saturated mantle peridotite to 11 gigapascals. Science 276, 240–243 (1997)

% Read the CSV file
filename = 'Kawamoto_Peridotite_water_saturated_solidus.csv';

% Read the data (skip the first row if it contains headers)
data = readmatrix(filename, 'NumHeaderLines', 1);

% Extract temperature (first column) and pressure (second column)
temperature = data(:, 1); % Temperature in °C
pressure = data(:, 2);    % Pressure in GPa

% Constants for depth calculation
rho = 3300; % Density of the upper mantle in kg/m^3 (approx.)
g = 9.81;   % Acceleration due to gravity in m/s^2

% Calculate depth from pressure
depth = (pressure * 1e9) / (rho * g) / 1000; % Depth in kilometers

% Create a high-quality plot
figure('Position', [100, 100, 900, 700], 'Color', 'white');

% Plot temperature vs pressure
plot(temperature, pressure, 'LineWidth', 3, 'Color', 'b', 'Marker', 'o', ...
     'MarkerSize', 4, 'MarkerFaceColor', 'blue');

% Set plot title and labels
xlabel('Temperature (°C)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Pressure (GPa)', 'FontSize', 14, 'FontWeight', 'bold');

% Set grid for better readability
grid on;
set(gca, 'GridAlpha', 0.3);

% Set axis limits to fit the data with some padding
xlim([600, 1400]);

% Reverse the direction of the y-axis for Pressure (increasing depth downward)
set(gca, 'YDir', 'reverse');

% Create secondary y-axis for Depth on the right side
yyaxis right;

% Plot the same data but with depth values (invisible line for axis scaling)
plot(temperature, depth, 'LineWidth', 0.1, 'Color', 'none');

% Set depth axis properties
ylabel('Depth (km)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'r');
set(gca, 'YColor', 'r'); % Red color for depth axis
set(gca, 'YDir', 'reverse'); % Reverse depth axis (increasing downward)

% Set appropriate depth axis limits
depth_min = min(depth);
depth_max = max(depth);
ylim([depth_min-0.1*(depth_max-depth_min), depth_max+0.1*(depth_max-depth_min)]);

% Switch back to left y-axis for future plotting
yyaxis left;

% Add legend
legend('Water-Saturated Peridotite Solidus (Kawamoto & Holloway, 1997)', 'Location', 'best', 'FontSize', 12);

% Display basic statistics in command window
fprintf('Kawamoto Water-Saturated Peridotite Solidus Data Summary:\n');
fprintf('Number of data points: %d\n', length(temperature));
fprintf('Temperature range: %.1f to %.1f °C\n', min(temperature), max(temperature));
fprintf('Pressure range: %.2f to %.2f GPa\n', min(pressure), max(pressure));
fprintf('Depth range: %.1f to %.1f km\n', min(depth), max(depth));

% Export the plot
saveas(gcf, 'kawamoto_water_saturated_solidus.png');
print('kawamoto_water_saturated_solidus_highres', '-dpng', '-r300');
savefig('kawamoto_water_saturated_solidus.fig');

fprintf('Plot saved successfully!\n');
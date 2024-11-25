% Define temperature range (1400 to 2400 K)
Temperature = 1500:10:2400; % in Kelvin
Pressure = 200; % Example pressure in kbar (constant for this plot)

% Initialize arrays for water solubilities and partition ratios
wad_solubility = zeros(size(Temperature));
ring_solubility = zeros(size(Temperature));
partition_wad_vs_ring = zeros(size(Temperature));
partition_wad_vs_gt = zeros(size(Temperature));

% Compute solubilities and partition ratios for each temperature
for i = 1:length(Temperature)
    % Call the water partition function
    output = water_partition_wad_ring_gt(Temperature(i) - 273.15, Pressure); % Convert K to °C
    wad_solubility(i) = output(1);        % Wadsleyite water solubility
    ring_solubility(i) = output(2);       % Ringwoodite water solubility
    partition_wad_vs_ring(i) = output(3); % Wad/Ring partition ratio
    partition_wad_vs_gt(i) = output(4);   % Wad/Garnet partition ratio
end


% Plot the data
figure('Color', 'white', 'Units', 'inches', 'Position', [1, 1, 10, 6]);

% Create first y-axis (water solubilities)
yyaxis left;
plot(Temperature, wad_solubility, '-o', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Wadsleyite Solubility');
hold on;
plot(Temperature, ring_solubility, '-s', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Ringwoodite Solubility');
ylabel('Water Solubility (wt%)', 'FontSize', 14, 'FontWeight', 'bold');
ylim([0, max([wad_solubility, ring_solubility]) * 1.1]); % Adjust y-axis range
set(gca, 'YColor', [0, 0.447, 0.741]); % Set left y-axis color
xlim([1500, 2400]);
ylim([0, 4]);
% Create second y-axis (partition ratios)
yyaxis right;
plot(Temperature, partition_wad_vs_ring, '-d', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Wadsleyite/Ringwoodite Partition');
hold on;
plot(Temperature, partition_wad_vs_gt, '-^', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Wadsleyite/Garnet Partition');
ylabel('Water Partition Ratios', 'FontSize', 14, 'FontWeight', 'bold');
ylim([0, max([partition_wad_vs_ring, partition_wad_vs_gt]) * 1.1]); % Adjust y-axis range
set(gca, 'YColor', [0.850, 0.325, 0.098]); % Set right y-axis color

% Add title and x-axis label
xlabel('Temperature (K)', 'FontSize', 14, 'FontWeight', 'bold');
title('Water Solubility and Partition Ratios vs Temperature', 'FontSize', 16, 'FontWeight', 'bold');

% General plot settings
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
legend('Location', 'best', 'FontSize', 12);

% Save the figure as a high-quality PNG file
exportgraphics(gcf, 'water_partition_wad_ring_gt_plot.png', 'Resolution', 300);

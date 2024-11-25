% Input Parameters
group_id = "KARATO";
uppermantle_layers = 29;
transiztionzone_layers = 27;
lowermantle_layers = 0;
fileName = "karato2011pyrolite_2_without_header_nocpxopxgrt.tab";

% Perform Calculations
[columnNames, dataMatrix, ~] = read_tab_format_without_header(fileName);
[upperMantleTable, ~, ~, ~] = divideDataByDominantMineral(dataMatrix, columnNames);
m = size(upperMantleTable, 1); % Number of rows
water_content = ones(m, 1);    % Water content (wt%)

[~, ~, sigma_plus_1, sigma_minus_1] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content, group_id);
[~, ~, sigma_plus_01, sigma_minus_01] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content * 0.1, group_id);
[~, ~, sigma_plus_001, sigma_minus_001] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content * 0.01, group_id);
[Pressure, Tempature, sigma_plus_0001, sigma_minus_0001] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content * 0.001, group_id);

% Convert Pressure to Depth (km) using a standard loop
depth = zeros(length(Pressure), 1); % Preallocate for efficiency
for i = 1:length(Pressure)
    depth(i) = find_depth_for_pressure(Pressure(i) * 0.0001); % Convert bar to GPa
end

% Validate Dimensions
assert(length(Pressure) == length(sigma_plus_1) && ...
       length(Pressure) == length(sigma_minus_1), ...
       'Error: All input vectors must have the same length.');

% Plot Temperature Profile
figure('Color', 'white', 'Units', 'inches', 'Position', [1, 1, 5, 6]);
plot(depth, Tempature, '-', 'LineWidth', 2, 'DisplayName', 'Temperature');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
ylim([1400, 2000]);
xlim([100, 400]);
xlabel('Depth (km)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('T (K)', 'FontSize', 14, 'FontWeight', 'bold');
title('Temperature Profile', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
box on;

% Plot Upper Bound Conductivities
figure('Color', 'white', 'Units', 'inches', 'Position', [1, 1, 5, 6]);
hold on;
plot(depth, sigma_plus_1, '-', 'LineWidth', 2, 'DisplayName', '1 wt%');
plot(depth, sigma_plus_01, '-', 'LineWidth', 2, 'DisplayName', '0.1 wt%');
plot(depth, sigma_plus_001, '-', 'LineWidth', 2, 'DisplayName', '0.01 wt%');
plot(depth, sigma_plus_0001, '-', 'LineWidth', 2, 'DisplayName', '0.001 wt%');
hold off;
set(gca, 'YScale', 'log', 'FontSize', 12, 'LineWidth', 1.5);
ylim([1e-4, 10]);
xlim([100, 400]);
xlabel('Depth (km)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Conductivity (S/m)', 'FontSize', 14, 'FontWeight', 'bold');
title('Upper Bound Conductivity', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
box on;

% Plot Lower Bound Conductivities
figure('Color', 'white', 'Units', 'inches', 'Position', [1, 1, 5, 6]);
hold on;
plot(depth, sigma_minus_1, '-', 'LineWidth', 2, 'DisplayName', '1 wt%');
plot(depth, sigma_minus_01, '-', 'LineWidth', 2, 'DisplayName', '0.1 wt%');
plot(depth, sigma_minus_001, '-', 'LineWidth', 2, 'DisplayName', '0.01 wt%');
plot(depth, sigma_minus_0001, '-', 'LineWidth', 2, 'DisplayName', '0.001 wt%');
hold off;
set(gca, 'YScale', 'log', 'FontSize', 12, 'LineWidth', 1.5);
ylim([1e-4, 10]);
xlim([100, 400]);
xlabel('Depth (km)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Conductivity (S/m)', 'FontSize', 14, 'FontWeight', 'bold');
title('Lower Bound Conductivity', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
box on;

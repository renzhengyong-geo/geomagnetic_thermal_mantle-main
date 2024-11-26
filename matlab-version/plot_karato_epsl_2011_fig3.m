% Input Parameters
group_id = "KARATO";
fileName = "../perple_x/karato2011pyrolite_2_without_header_nocpxopxgrt.tab"; % 3.1Gpa to 23.0 Gpa
% fileName = "../perple_x/karato2011pyrolite_2.tab"; % 3.1Gpa to 23.0 Gpa

% Perform Calculations
[columnNames, dataMatrix, ~] = read_tab_format_without_header(fileName);
[upperMantleTable, upperTransitionZoneTable, lowerTransitionZoneTable, ~] = divideDataByDominantMineral(dataMatrix, columnNames);

% Compute conductivity-depth profile in the upper mantle
m = size(upperMantleTable, 1); % Number of rows
water_content = ones(m, 1);    % Water content (wt%)
[~, ~, um_sigma_plus_1, um_sigma_minus_1, um_Cw1] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content, group_id);
[~, ~, um_sigma_plus_01, um_sigma_minus_01, um_Cw01] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content * 0.1, group_id);
[~, ~, um_sigma_plus_001, um_sigma_minus_001, um_Cw001] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content * 0.01, group_id);
[um_P, um_T, um_sigma_plus_0001, um_sigma_minus_0001,um_Cw0001] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content * 0.001, group_id);

% Compute conductivity-depth profile in the transition zone
m = size(upperTransitionZoneTable, 1); % Number of rows
water_content = ones(m, 1);    % Water content (wt%)
[~, ~, ut_sigma_plus_1, ut_sigma_minus_1, ut_Cw1] = compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content, group_id);
[~, ~, ut_sigma_plus_01, ut_sigma_minus_01,ut_Cw01] = compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content * 0.1, group_id);
[~, ~, ut_sigma_plus_001, ut_sigma_minus_001,ut_Cw001] = compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content * 0.01, group_id);
[ut_P, ut_T, ut_sigma_plus_0001, ut_sigma_minus_0001, ut_Cw0001] = compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content * 0.001, group_id);

m = size(lowerTransitionZoneTable, 1); % Number of rows
water_content = ones(m, 1);    % Water content (wt%)
[~, ~, lt_sigma_plus_1, lt_sigma_minus_1, lt_Cw1] = compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content, group_id);
[~, ~, lt_sigma_plus_01, lt_sigma_minus_01,lt_Cw01] = compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content * 0.1, group_id);
[~, ~, lt_sigma_plus_001, lt_sigma_minus_001,lt_Cw001] = compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content * 0.01, group_id);
[lt_P, lt_T, lt_sigma_plus_0001, lt_sigma_minus_0001, lt_Cw0001] = compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content * 0.001, group_id);

% Assemble P, T, and sigma
P = [um_P; ut_P; lt_P];
T = [um_T; ut_T; lt_T];
Cw1=[um_Cw1;ut_Cw1;lt_Cw1];
Cw01=[um_Cw01;ut_Cw01;lt_Cw01];
Cw001=[um_Cw001;ut_Cw001;lt_Cw001];
Cw0001=[um_Cw0001;ut_Cw0001;lt_Cw0001];
% Display water content arrays for verification
disp('Cw1 (1 wt% water content):');
disp(Cw1);
disp('Cw01 (0.1 wt% water content):');
disp(Cw01);
disp('Cw001 (0.01 wt% water content):');
disp(Cw001);
disp('Cw0001 (0.001 wt% water content):');
disp(Cw0001);

sigma_plus_1 = [um_sigma_plus_1; ut_sigma_plus_1; lt_sigma_plus_1];
sigma_plus_01 = [um_sigma_plus_01; ut_sigma_plus_01; lt_sigma_plus_01];
sigma_plus_001 = [um_sigma_plus_001; ut_sigma_plus_001; lt_sigma_plus_001];
sigma_plus_0001 = [um_sigma_plus_0001; ut_sigma_plus_0001; lt_sigma_plus_0001];
sigma_minus_1 = [um_sigma_minus_1; ut_sigma_minus_1; lt_sigma_minus_1];
sigma_minus_01 = [um_sigma_minus_01; ut_sigma_minus_01; lt_sigma_minus_01];
sigma_minus_001 = [um_sigma_minus_001; ut_sigma_minus_001; lt_sigma_minus_001];
sigma_minus_0001 = [um_sigma_minus_0001; ut_sigma_minus_0001; lt_sigma_minus_0001];

% Convert Pressure to Depth (km) using a standard loop
depth = zeros(length(P), 1); % Preallocate for efficiency
for i = 1:length(P)
    depth(i) = find_depth_for_pressure(P(i) * 0.0001); % Convert bar to GPa
end

% Plot Temperature Profile
figure('Color', 'white', 'Units', 'inches', 'Position', [1, 1, 5, 6]);
plot(depth, T, '-', 'LineWidth', 2, 'DisplayName', 'Temperature');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
ylim([1400, 2000]);
xlim([100, 660]);
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
xlim([100, 660]);
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
xlim([100, 660]);
xlabel('Depth (km)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Conductivity (S/m)', 'FontSize', 14, 'FontWeight', 'bold');
title('Lower Bound Conductivity', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
box on;

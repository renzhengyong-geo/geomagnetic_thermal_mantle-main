
% Constants
T = 1573;              % Temperature in Kelvin
P = 2e9;               % Pressure in Pa
Xfe = 1.0;             % Fe content (dimensionless)
Cw_values = logspace(-6, 2, 100); % Water content range in wt% (log scale from 0.0001 to 1)

% Initialize arrays to store conductivity values
sigma_Yoshino = zeros(size(Cw_values));
sigma_Karato = zeros(size(Cw_values));

% Calculate conductivity for each water content value
for i = 1:length(Cw_values)
    Cw = Cw_values(i);
    % Use Yoshino olivine conductivity subroutine
    sigma_Yoshino(i) = Yoshino_olivine_conductivity(T, Cw, P, Xfe);
    % Use Karato olivine conductivity subroutine
    sigma_Karato(i) = Karato_olivine_conductivity(T, Cw, P);
end

% Plot the results with log-log scale
figure('Units', 'inches', 'Position', [0 0 6 6.5 ], 'PaperPositionMode', 'auto');

loglog(Cw_values, sigma_Yoshino, '--r', 'LineWidth', 4); % Yoshino curve
hold on;
loglog(Cw_values, sigma_Karato, '--k', 'LineWidth', 4); % Karato curve
hold off;

% Add labels, legend, and title
xlabel('Water Content (wt%)');
ylabel('Conductivity (S/m)');
title('Olivine Conductivity via Water Content');
legend('Yoshino Olivine Conductivity', 'Wang and Karato Olivine Conductivity', 'Location', 'Best');
grid on;

% Set axis limits and ticks for better visualization
xlim([1e-4, 1e0]); % Water content from 0.01 wt% to 100 wt%
ylim([1e-3, 1e2]); % Conductivity range (adjust based on expected results)
set(gca, 'FontSize', 14); % Increase font size for better readability

% Save high-quality plot
set(gcf, 'PaperPositionMode', 'auto'); % Use screen size for printing
print('Yosino_Olivine_Conductivity', '-dpng', '-r300'); % Save as PNG with 300 dpi


% Clear workspace and initialize
clear;
clc;
% Constants for both plots
P = 15e9;                   % Pressure in Pa (15 GPa)
Xfe = 0.1;                  % Fe content for Yoshino model (dimensionless)
T_values = linspace(500, 2000, 100); % Temperature range in K
Cw_values = linspace(1e-3, 5, 100);   % Water content range in wt% for second plot

% Initialize arrays for conductivity values
sigma_Yoshino_0wt = zeros(size(T_values));
sigma_Yoshino_0_1wt = zeros(size(T_values));
sigma_Karato_0wt = zeros(size(T_values));
sigma_Karato_0_1wt = zeros(size(T_values));

% Plot 1: Conductivity vs. 1000/T
for i = 1:length(T_values)
    T = T_values(i); % Use T directly in Kelvin
    sigma_Yoshino_0wt(i) = Yoshino_wadsleyite_conductivity(T, 0, P, Xfe);    % Cw = 0 wt%
    sigma_Yoshino_0_1wt(i) = Yoshino_wadsleyite_conductivity(T, 0.1, P, Xfe); % Cw = 0.1 wt%
    sigma_Karato_0wt(i) = Karato_wadsleyite_conductivity(T, 0, P);          % Cw = 0 wt%
    sigma_Karato_0_1wt(i) = Karato_wadsleyite_conductivity(T, 0.1, P);      % Cw = 0.1 wt%
end

figure('Units', 'inches', 'Position', [0 0 6 6.5 ], 'PaperPositionMode', 'auto');
% Define RGB color for #EDB120
color_EDB120 = [237, 177, 32] / 255; % Convert from 0-255 range to 0-1

% Plot with the custom color
semilogy(1000 ./ T_values, sigma_Yoshino_0wt, '-r', 'LineWidth', 2); hold on;
semilogy(1000 ./ T_values, sigma_Yoshino_0_1wt, '--r', 'LineWidth', 2);
semilogy(1000 ./ T_values, sigma_Karato_0wt, '-', 'Color', color_EDB120, 'LineWidth', 2);
semilogy(1000 ./ T_values, sigma_Karato_0_1wt, '--', 'Color', color_EDB120, 'LineWidth', 2);
hold off;

% Customize the first plot
xlabel('1000 / T (K^{-1})');
ylabel('Conductivity (\sigma) [S/m]');
title('Wadsleyite Conductivity vs. 1000/T');
legend('Yoshino (Cw = 0 wt%)', 'Yoshino (Cw = 0.1 wt%)', ...
       'Karato (Cw = 0 wt%)', 'Karato (Cw = 0.1 wt%)', 'Location', 'Best');
grid on;
xlim([0.5, 1.0]); % Adjusted for 1000/T range
ylim([1e-4, 1e0]); % Conductivity range
set(gca, 'FontSize', 14); % Increase font size for readability

% Save first plot as high-quality image
set(gcf, 'PaperPositionMode', 'auto');
print('Yosino_wadsleyite_conductivity_vs_1000_T_4Curves', '-dpng', '-r300');

% Plot 2: Conductivity vs. Water Content
T_fixed = 1700; % Fixed temperature for second plot
sigma_Yoshino_Cw = zeros(size(Cw_values));
sigma_Karato_Cw = zeros(size(Cw_values));

for i = 1:length(Cw_values)
    sigma_Yoshino_Cw(i) = Yoshino_wadsleyite_conductivity(T_fixed, Cw_values(i), P, Xfe);
    sigma_Karato_Cw(i) = Karato_wadsleyite_conductivity(T_fixed, Cw_values(i), P);
end

% Create second plot: Conductivity vs. Water Content
figure('Units', 'inches', 'Position', [0 0 6 6.5 ], 'PaperPositionMode', 'auto');
loglog(Cw_values, sigma_Yoshino_Cw, '-r', 'LineWidth', 2); hold on;
loglog(Cw_values, sigma_Karato_Cw, '-', 'Color', color_EDB120, 'LineWidth', 2);
hold off;

% Customize the second plot
xlabel('Water Content (wt%)');
ylabel('Conductivity (\sigma) [S/m]');
title('Wadsleyite at 1700 K');
legend('Yoshino', 'Karato', 'Location', 'Best');
grid on;
xlim([1e-3, 2]);
ylim([1e-4, 10]);
set(gca, 'FontSize', 14); % Increase font size for readability

% Save second plot as high-quality image
set(gcf, 'PaperPositionMode', 'auto');
print('Yosino_wadsleyite_conductivity_vs_WaterContent', '-dpng', '-r300');


% Clear workspace and initialize
clear;
clc;
% Constants for both plots
P = 15e9;                   % Pressure in Pa (15 GPa)
Xfe = 1.0;                  % Fe content for Yoshino model (dimensionless)
T_values = linspace(500, 2000, 100); % Temperature range in K
Cw_values = linspace(1e-3, 5, 100);   % Water content range in wt% for second plot

% Initialize arrays for conductivity values
sigma_Yoshino_0wt = zeros(size(T_values));
sigma_Yoshino_0_1wt = zeros(size(T_values));
sigma_Karato_0wt = zeros(size(T_values));
sigma_Karato_0_1wt = zeros(size(T_values));

% Plot 1: Conductivity vs. 1000/T
for i = 1:length(T_values)
    T = T_values(i); % Use T directly in Kelvin
    sigma_Yoshino_0wt(i) = Yoshino_ringwoodite_conductivity(T, 0, P, Xfe);    % Cw = 0 wt%
    sigma_Yoshino_0_1wt(i) = Yoshino_ringwoodite_conductivity(T, 0.1, P, Xfe); % Cw = 0.1 wt%   
    sigma_Karato_0wt(i) = Karato_ringwoodite_conductivity(T, 0.001, P);          % Cw = 0 wt%
    sigma_Karato_0_1wt(i) = Karato_ringwoodite_conductivity(T, 0.1, P);      % Cw = 0.1 wt%
end

figure('Units', 'inches', 'Position', [0 0 6 6.5 ], 'PaperPositionMode', 'auto');
% Define RGB color for #EDB120
color_EDB120 = [237, 177, 32] / 255; % Convert from 0-255 range to 0-1

% Plot with the custom color
semilogy(1000 ./ T_values, sigma_Yoshino_0wt, '-r', 'LineWidth', 2); hold on;
semilogy(1000 ./ T_values, sigma_Yoshino_0_1wt, '--r', 'LineWidth', 2);
semilogy(1000 ./ T_values, sigma_Karato_0wt, '-', 'Color', color_EDB120, 'LineWidth', 2);
semilogy(1000 ./ T_values, sigma_Karato_0_1wt, '--', 'Color', color_EDB120, 'LineWidth', 2);
hold off;

% Customize the first plot
xlabel('1000 / T (K^{-1})');
ylabel('Conductivity (\sigma) [S/m]');
title('ringwoodite Conductivity vs. 1000/T');
legend('Yoshino (Cw = 0 wt%)', 'Yoshino (Cw = 0.1 wt%)', ...
       'Karato (Cw = 0 wt%)', 'Karato (Cw = 0.1 wt%)', 'Location', 'Best');
grid on;
xlim([0.5, 1.0]); % Adjusted for 1000/T range
ylim([1e-4, 1e0]); % Conductivity range
set(gca, 'FontSize', 14); % Increase font size for readability

% Save first plot as high-quality image
set(gcf, 'PaperPositionMode', 'auto');
print('Yosino_ringwoodite_conductivity_vs_1000_T_4Curves', '-dpng', '-r300');

% Plot 2: Conductivity vs. Water Content
T_fixed = 1700; % Fixed temperature for second plot
sigma_Yoshino_Cw = zeros(size(Cw_values));
sigma_Karato_Cw = zeros(size(Cw_values));

for i = 1:length(Cw_values)
    sigma_Yoshino_Cw(i) = Yoshino_ringwoodite_conductivity(T_fixed, Cw_values(i), P, Xfe);
    sigma_Karato_Cw(i) = Karato_ringwoodite_conductivity(T_fixed, Cw_values(i), P);
end

% Create second plot: Conductivity vs. Water Content
figure('Units', 'inches', 'Position', [0 0 6 6.5 ], 'PaperPositionMode', 'auto');
loglog(Cw_values, sigma_Yoshino_Cw, '-r', 'LineWidth', 2); hold on;
loglog(Cw_values, sigma_Karato_Cw, '-', 'Color', color_EDB120, 'LineWidth', 2);
hold off;

% Customize the second plot
xlabel('Water Content (wt%)');
ylabel('Conductivity (\sigma) [S/m]');
title('ringwoodite at 1700 K');
legend('Yoshino', 'Karato', 'Location', 'Best');
grid on;
xlim([1e-3, 2]);
ylim([1e-4, 10]);
set(gca, 'FontSize', 14); % Increase font size for readability

% Save second plot as high-quality image
set(gcf, 'PaperPositionMode', 'auto');
print('Yosino_ringwoodite_Conductivity_vs_WaterContent', '-dpng', '-r300');

% Plot the condudviity-depth profile using the average mantle geotherm
% (Katsura et al.2010, PEPI）
% Clear workspace and initialize
clear;
clc;
% % Define the filename
% filename = 'geothermal_katsura_PEPI_2010.csv'; % Replace with your file name
% % Read the CSV file, skipping the first row
% opts = detectImportOptions(filename);  % Automatically detect the file structure
% opts.DataLines = [2, Inf];  % Start reading from the second row
% data = readtable(filename, opts);  % Read the table
% % Extract depth and temperature columns
% depth = data.Depth;  % Adjust if your column names are different
% temperature = data.T;  % Adjust if your column names are different
% pressure=zeros(length(depth),1); 
% for i=1:length(depth)
%     pressure(i)=find_pressure_for_depth(depth(i)); %% input in unit of km
% end
% pressure=pressure*1e4; %%GPa to Bar
% disp([pressure', temperature', depth]);
% 
% Input Parameters
group_id = "YoSHINO";
fileName = "../perple_x/Yoshino_AREPS_2013_pyrolite.tab"; % for fig 7

% Perform Calculations
[columnNames, dataMatrix, ~] = read_tab_format_without_header(fileName);
[upperMantleTable, upperTransitionZoneTable, lowerTransitionZoneTable, ~] = divideDataByDominantMineral(dataMatrix, columnNames);
% Compute conductivity-depth profile in the upper mantle
m = size(upperMantleTable, 1); % Number of rows
water_content = ones(m, 1);    % Water content (wt%)
[~, ~, um_sigma_plus_1, um_sigma_minus_1, um_Cw1] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content, group_id);
[~, ~, um_sigma_plus_01, um_sigma_minus_01, um_Cw01] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content * 0.5, group_id);
[~, ~, um_sigma_plus_001, um_sigma_minus_001, um_Cw001] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content * 0.1, group_id);
[um_P, um_T, um_sigma_plus_0001, um_sigma_minus_0001,um_Cw0001] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content * 1e-5, group_id);
% Compute conductivity-depth profile in the transition zone
m = size(upperTransitionZoneTable, 1); % Number of rows
water_content = ones(m, 1);    % Water content (wt%)
[~, ~, ut_sigma_plus_1, ut_sigma_minus_1, ut_Cw1] = compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content, group_id);
[~, ~, ut_sigma_plus_01, ut_sigma_minus_01,ut_Cw01] = compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content * 0.5, group_id);
[~, ~, ut_sigma_plus_001, ut_sigma_minus_001,ut_Cw001] = compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content * 0.1, group_id);
[ut_P, ut_T, ut_sigma_plus_0001, ut_sigma_minus_0001, ut_Cw0001] = compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content * 1e-5, group_id);
m = size(lowerTransitionZoneTable, 1); % Number of rows
water_content = ones(m, 1);    % Water content (wt%)
[~, ~, lt_sigma_plus_1, lt_sigma_minus_1, lt_Cw1] = compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content, group_id);
[~, ~, lt_sigma_plus_01, lt_sigma_minus_01,lt_Cw01] = compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content * 0.5, group_id);
[~, ~, lt_sigma_plus_001, lt_sigma_minus_001,lt_Cw001] = compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content * 0.1, group_id);
[lt_P, lt_T, lt_sigma_plus_0001, lt_sigma_minus_0001, lt_Cw0001] = compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content * 1e-5, group_id);

% Assemble P, T, and sigma
P = [um_P; ut_P; lt_P];
T = [um_T; ut_T; lt_T];
Cw1=[um_Cw1;ut_Cw1;lt_Cw1];
Cw01=[um_Cw01;ut_Cw01;lt_Cw01];
Cw001=[um_Cw001;ut_Cw001;lt_Cw001];
Cw0001=[um_Cw0001;ut_Cw0001;lt_Cw0001];
% Display water content arrays for verification
disp('Cw1 (1.0 wt% water content):');
disp(Cw1);
disp('Cw01 (0.5 wt% water content):');
disp(Cw01);
disp('Cw001 (0.1 wt% water content):');
disp(Cw001);
disp('Cw0001 (0.0 wt% water content):');
disp(Cw0001);

sigma_plus_1 = [um_sigma_plus_1; ut_sigma_plus_1; lt_sigma_plus_1];
sigma_plus_01 = [um_sigma_plus_01; ut_sigma_plus_01; lt_sigma_plus_01];
sigma_plus_001 = [um_sigma_plus_001; ut_sigma_plus_001; lt_sigma_plus_001];
sigma_plus_0001 = [um_sigma_plus_0001; ut_sigma_plus_0001; lt_sigma_plus_0001];
sigma_minus_1 = [um_sigma_minus_1; ut_sigma_minus_1; lt_sigma_minus_1];
sigma_minus_01 = [um_sigma_minus_01; ut_sigma_minus_01; lt_sigma_minus_01];
sigma_minus_001 = [um_sigma_minus_001; ut_sigma_minus_001; lt_sigma_minus_001];
sigma_minus_0001 = [um_sigma_minus_0001; ut_sigma_minus_0001; lt_sigma_minus_0001];

y_sigma_plus_1 =sigma_plus_1;
y_sigma_plus_01=sigma_plus_01;
y_sigma_plus_001=sigma_plus_001;
y_sigma_plus_0001=sigma_plus_0001;
y_sigma_minus_1=sigma_minus_1;
y_sigma_minus_01=sigma_minus_01;
y_sigma_minus_001=sigma_minus_001;
y_sigma_minus_0001=sigma_minus_0001;


group_id = "KARATO";
fileName = "../perple_x/Yoshino_AREPS_2013_pyrolite.tab"; % for fig 7

% Perform Calculations
[columnNames, dataMatrix, ~] = read_tab_format_without_header(fileName);
[upperMantleTable, upperTransitionZoneTable, lowerTransitionZoneTable, ~] = divideDataByDominantMineral(dataMatrix, columnNames);
% Compute conductivity-depth profile in the upper mantle
m = size(upperMantleTable, 1); % Number of rows
water_content = ones(m, 1);    % Water content (wt%)
[~, ~, um_sigma_plus_1, um_sigma_minus_1, um_Cw1] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content, group_id);
[~, ~, um_sigma_plus_01, um_sigma_minus_01, um_Cw01] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content * 0.1, group_id);
[~, ~, um_sigma_plus_001, um_sigma_minus_001, um_Cw001] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content * 0.01, group_id);
[um_P, um_T, um_sigma_plus_0001, um_sigma_minus_0001,um_Cw0001] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content * 1e-5, group_id);
% Compute conductivity-depth profile in the transition zone
m = size(upperTransitionZoneTable, 1); % Number of rows
water_content = ones(m, 1);    % Water content (wt%)
[~, ~, ut_sigma_plus_1, ut_sigma_minus_1, ut_Cw1] = compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content, group_id);
[~, ~, ut_sigma_plus_01, ut_sigma_minus_01,ut_Cw01] = compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content * 0.1, group_id);
[~, ~, ut_sigma_plus_001, ut_sigma_minus_001,ut_Cw001] = compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content * 0.01, group_id);
[ut_P, ut_T, ut_sigma_plus_0001, ut_sigma_minus_0001, ut_Cw0001] = compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content * 1e-5, group_id);
m = size(lowerTransitionZoneTable, 1); % Number of rows
water_content = ones(m, 1);    % Water content (wt%)
[~, ~, lt_sigma_plus_1, lt_sigma_minus_1, lt_Cw1] = compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content, group_id);
[~, ~, lt_sigma_plus_01, lt_sigma_minus_01,lt_Cw01] = compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content * 0.1, group_id);
[~, ~, lt_sigma_plus_001, lt_sigma_minus_001,lt_Cw001] = compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content * 0.01, group_id);
[lt_P, lt_T, lt_sigma_plus_0001, lt_sigma_minus_0001, lt_Cw0001] = compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content * 1e-5, group_id);

sigma_plus_1 = [um_sigma_plus_1; ut_sigma_plus_1; lt_sigma_plus_1];
sigma_plus_01 = [um_sigma_plus_01; ut_sigma_plus_01; lt_sigma_plus_01];
sigma_plus_001 = [um_sigma_plus_001; ut_sigma_plus_001; lt_sigma_plus_001];
sigma_plus_0001 = [um_sigma_plus_0001; ut_sigma_plus_0001; lt_sigma_plus_0001];
sigma_minus_1 = [um_sigma_minus_1; ut_sigma_minus_1; lt_sigma_minus_1];
sigma_minus_01 = [um_sigma_minus_01; ut_sigma_minus_01; lt_sigma_minus_01];
sigma_minus_001 = [um_sigma_minus_001; ut_sigma_minus_001; lt_sigma_minus_001];
sigma_minus_0001 = [um_sigma_minus_0001; ut_sigma_minus_0001; lt_sigma_minus_0001];

k_sigma_plus_1 =sigma_plus_1;
k_sigma_plus_01=sigma_plus_01;
k_sigma_plus_001=sigma_plus_001;
k_sigma_plus_0001=sigma_plus_0001;
k_sigma_minus_1=sigma_minus_1;
k_sigma_minus_01=sigma_minus_01;
k_sigma_minus_001=sigma_minus_001;
k_sigma_minus_0001=sigma_minus_0001;

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
xlim([200, 700]);
xlabel('Depth (km)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('T (K)', 'FontSize', 14, 'FontWeight', 'bold');
title('Temperature Profile', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
box on;

% Plot Upper Bound Conductivities
figure('Color', 'white', 'Units', 'inches', 'Position', [1, 1, 6, 11]);
subplot(2,1,1);
plot(depth, y_sigma_plus_1, '-', 'LineWidth', 2, 'DisplayName', '1 wt%');
hold on;
plot(depth, y_sigma_plus_01, '-', 'LineWidth', 2, 'DisplayName', '0.5 wt%');
plot(depth, y_sigma_plus_001, '-', 'LineWidth', 2, 'DisplayName', '0.1 wt%');
plot(depth, y_sigma_plus_0001, '-', 'LineWidth', 2, 'DisplayName', '1e-5 wt%');
hold off;
set(gca, 'YScale', 'log', 'FontSize', 12, 'LineWidth', 1.5);
ylim([1e-3, 10]);
xlim([200, 700]);
% xlabel('Depth (km)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Conductivity (S/m)', 'FontSize', 14, 'FontWeight', 'bold');
title('Yoshino model - HS upper bound', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
box on;
subplot(2,1,2);
plot(depth, k_sigma_plus_1, '-', 'LineWidth', 2, 'DisplayName', '1 wt%');
hold on;
plot(depth, k_sigma_plus_01, '-', 'LineWidth', 2, 'DisplayName', '0.1 wt%');
plot(depth, k_sigma_plus_001, '-', 'LineWidth', 2, 'DisplayName', '0.01 wt%');
plot(depth, k_sigma_plus_0001, '-', 'LineWidth', 2, 'DisplayName', '1e-5 wt%');
hold off;
set(gca, 'YScale', 'log', 'FontSize', 12, 'LineWidth', 1.5);
ylim([1e-3, 10]);
xlim([200, 700]);
xlabel('Depth (km)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Conductivity (S/m)', 'FontSize', 14, 'FontWeight', 'bold');
title('Karato model - HS upper bound', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
box on;

% Export the plot to a high-quality PNG file
exportgraphics(gcf, 'Yoshino_AREPS_2013_fig7_ub.png', 'Resolution', 300);

% Plot Lower Bound Conductivities
figure('Color', 'white', 'Units', 'inches', 'Position', [1, 1, 6, 11]);
subplot(2,1,1);
plot(depth, y_sigma_minus_1, '-', 'LineWidth', 2, 'DisplayName', '1 wt%');
hold on;
plot(depth, y_sigma_minus_01, '-', 'LineWidth', 2, 'DisplayName', '0.5 wt%');
plot(depth, y_sigma_minus_001, '-', 'LineWidth', 2, 'DisplayName', '0.1 wt%');
plot(depth, y_sigma_minus_0001, '-', 'LineWidth', 2, 'DisplayName', '1e-5 wt%');
hold off;
set(gca, 'YScale', 'log', 'FontSize', 12, 'LineWidth', 1.5);
ylim([1e-3, 10]);
xlim([200, 700]);
% xlabel('Depth (km)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Conductivity (S/m)', 'FontSize', 14, 'FontWeight', 'bold');
title('Lower Bound Conductivity', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
box on;
subplot(2,1,2);
plot(depth, k_sigma_minus_1, '-', 'LineWidth', 2, 'DisplayName', '1 wt%');
hold on;
plot(depth, k_sigma_minus_01, '-', 'LineWidth', 2, 'DisplayName', '0.1 wt%');
plot(depth, k_sigma_minus_001, '-', 'LineWidth', 2, 'DisplayName', '0.01 wt%');
plot(depth, k_sigma_minus_0001, '-', 'LineWidth', 2, 'DisplayName', '1e-5 wt%');
hold off;
set(gca, 'YScale', 'log', 'FontSize', 12, 'LineWidth', 1.5);
ylim([1e-3, 10]);
xlim([200, 700]);
xlabel('Depth (km)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Conductivity (S/m)', 'FontSize', 14, 'FontWeight', 'bold');
title('Lower Bound Conductivity', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
box on;
% Export the plot to a high-quality PNG file
exportgraphics(gcf, 'Yoshino_AREPS_2013_fig7_lb.png', 'Resolution', 300);


% Plot the condudviity-depth profile using the average mantle geotherm
% (Katsura et al.2010, PEPI）300k lower
% Clear workspace and initialize
clear;
clc;
% Input Parameters
group_id = "YoSHINO";
fileName = "../perple_x/Yoshino_AREPS_2013_pyrolite_300k_lower.tab"; % for fig 7

% Perform Calculations
[columnNames, dataMatrix, ~] = read_tab_format_without_header(fileName);
[upperMantleTable, upperTransitionZoneTable, lowerTransitionZoneTable, ~] = divideDataByDominantMineral(dataMatrix, columnNames);
% Compute conductivity-depth profile in the upper mantle
m = size(upperMantleTable, 1); % Number of rows
water_content = ones(m, 1);    % Water content (wt%)
[~, ~, um_sigma_plus_1, um_sigma_minus_1, um_Cw1] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content, group_id);
[~, ~, um_sigma_plus_01, um_sigma_minus_01, um_Cw01] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content * 0.5, group_id);
[~, ~, um_sigma_plus_001, um_sigma_minus_001, um_Cw001] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content * 0.1, group_id);
[um_P, um_T, um_sigma_plus_0001, um_sigma_minus_0001,um_Cw0001] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content * 1e-5, group_id);
% Compute conductivity-depth profile in the transition zone
m = size(upperTransitionZoneTable, 1); % Number of rows
water_content = ones(m, 1);    % Water content (wt%)
[~, ~, ut_sigma_plus_1, ut_sigma_minus_1, ut_Cw1] = compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content, group_id);
[~, ~, ut_sigma_plus_01, ut_sigma_minus_01,ut_Cw01] = compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content * 0.5, group_id);
[~, ~, ut_sigma_plus_001, ut_sigma_minus_001,ut_Cw001] = compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content * 0.1, group_id);
[ut_P, ut_T, ut_sigma_plus_0001, ut_sigma_minus_0001, ut_Cw0001] = compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content * 1e-5, group_id);
m = size(lowerTransitionZoneTable, 1); % Number of rows
water_content = ones(m, 1);    % Water content (wt%)
[~, ~, lt_sigma_plus_1, lt_sigma_minus_1, lt_Cw1] = compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content, group_id);
[~, ~, lt_sigma_plus_01, lt_sigma_minus_01,lt_Cw01] = compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content * 0.5, group_id);
[~, ~, lt_sigma_plus_001, lt_sigma_minus_001,lt_Cw001] = compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content * 0.1, group_id);
[lt_P, lt_T, lt_sigma_plus_0001, lt_sigma_minus_0001, lt_Cw0001] = compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content * 1e-5, group_id);

% Assemble P, T, and sigma
P = [um_P; ut_P; lt_P];
T = [um_T; ut_T; lt_T];
Cw1=[um_Cw1;ut_Cw1;lt_Cw1];
Cw01=[um_Cw01;ut_Cw01;lt_Cw01];
Cw001=[um_Cw001;ut_Cw001;lt_Cw001];
Cw0001=[um_Cw0001;ut_Cw0001;lt_Cw0001];
% Display water content arrays for verification
disp('Cw1 (1.0 wt% water content):');
disp(Cw1);
disp('Cw01 (0.5 wt% water content):');
disp(Cw01);
disp('Cw001 (0.1 wt% water content):');
disp(Cw001);
disp('Cw0001 (0.0 wt% water content):');
disp(Cw0001);

sigma_plus_1 = [um_sigma_plus_1; ut_sigma_plus_1; lt_sigma_plus_1];
sigma_plus_01 = [um_sigma_plus_01; ut_sigma_plus_01; lt_sigma_plus_01];
sigma_plus_001 = [um_sigma_plus_001; ut_sigma_plus_001; lt_sigma_plus_001];
sigma_plus_0001 = [um_sigma_plus_0001; ut_sigma_plus_0001; lt_sigma_plus_0001];
sigma_minus_1 = [um_sigma_minus_1; ut_sigma_minus_1; lt_sigma_minus_1];
sigma_minus_01 = [um_sigma_minus_01; ut_sigma_minus_01; lt_sigma_minus_01];
sigma_minus_001 = [um_sigma_minus_001; ut_sigma_minus_001; lt_sigma_minus_001];
sigma_minus_0001 = [um_sigma_minus_0001; ut_sigma_minus_0001; lt_sigma_minus_0001];

y_sigma_plus_1 =sigma_plus_1;
y_sigma_plus_01=sigma_plus_01;
y_sigma_plus_001=sigma_plus_001;
y_sigma_plus_0001=sigma_plus_0001;
y_sigma_minus_1=sigma_minus_1;
y_sigma_minus_01=sigma_minus_01;
y_sigma_minus_001=sigma_minus_001;
y_sigma_minus_0001=sigma_minus_0001;


group_id = "KARATO";
fileName = "../perple_x/Yoshino_AREPS_2013_pyrolite_300k_lower.tab"; % for fig 7

% Perform Calculations
[columnNames, dataMatrix, ~] = read_tab_format_without_header(fileName);
[upperMantleTable, upperTransitionZoneTable, lowerTransitionZoneTable, ~] = divideDataByDominantMineral(dataMatrix, columnNames);
% Compute conductivity-depth profile in the upper mantle
m = size(upperMantleTable, 1); % Number of rows
water_content = ones(m, 1);    % Water content (wt%)
[~, ~, um_sigma_plus_1, um_sigma_minus_1, um_Cw1] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content, group_id);
[~, ~, um_sigma_plus_01, um_sigma_minus_01, um_Cw01] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content * 0.1, group_id);
[~, ~, um_sigma_plus_001, um_sigma_minus_001, um_Cw001] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content * 0.01, group_id);
[um_P, um_T, um_sigma_plus_0001, um_sigma_minus_0001,um_Cw0001] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, water_content * 1e-5, group_id);
% Compute conductivity-depth profile in the transition zone
m = size(upperTransitionZoneTable, 1); % Number of rows
water_content = ones(m, 1);    % Water content (wt%)
[~, ~, ut_sigma_plus_1, ut_sigma_minus_1, ut_Cw1] = compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content, group_id);
[~, ~, ut_sigma_plus_01, ut_sigma_minus_01,ut_Cw01] = compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content * 0.1, group_id);
[~, ~, ut_sigma_plus_001, ut_sigma_minus_001,ut_Cw001] = compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content * 0.01, group_id);
[ut_P, ut_T, ut_sigma_plus_0001, ut_sigma_minus_0001, ut_Cw0001] = compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content * 1e-5, group_id);
m = size(lowerTransitionZoneTable, 1); % Number of rows
water_content = ones(m, 1);    % Water content (wt%)
[~, ~, lt_sigma_plus_1, lt_sigma_minus_1, lt_Cw1] = compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content, group_id);
[~, ~, lt_sigma_plus_01, lt_sigma_minus_01,lt_Cw01] = compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content * 0.1, group_id);
[~, ~, lt_sigma_plus_001, lt_sigma_minus_001,lt_Cw001] = compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content * 0.01, group_id);
[lt_P, lt_T, lt_sigma_plus_0001, lt_sigma_minus_0001, lt_Cw0001] = compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content * 1e-5, group_id);

sigma_plus_1 = [um_sigma_plus_1; ut_sigma_plus_1; lt_sigma_plus_1];
sigma_plus_01 = [um_sigma_plus_01; ut_sigma_plus_01; lt_sigma_plus_01];
sigma_plus_001 = [um_sigma_plus_001; ut_sigma_plus_001; lt_sigma_plus_001];
sigma_plus_0001 = [um_sigma_plus_0001; ut_sigma_plus_0001; lt_sigma_plus_0001];
sigma_minus_1 = [um_sigma_minus_1; ut_sigma_minus_1; lt_sigma_minus_1];
sigma_minus_01 = [um_sigma_minus_01; ut_sigma_minus_01; lt_sigma_minus_01];
sigma_minus_001 = [um_sigma_minus_001; ut_sigma_minus_001; lt_sigma_minus_001];
sigma_minus_0001 = [um_sigma_minus_0001; ut_sigma_minus_0001; lt_sigma_minus_0001];

k_sigma_plus_1 =sigma_plus_1;
k_sigma_plus_01=sigma_plus_01;
k_sigma_plus_001=sigma_plus_001;
k_sigma_plus_0001=sigma_plus_0001;
k_sigma_minus_1=sigma_minus_1;
k_sigma_minus_01=sigma_minus_01;
k_sigma_minus_001=sigma_minus_001;
k_sigma_minus_0001=sigma_minus_0001;

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
xlim([200, 700]);
xlabel('Depth (km)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('T (K)', 'FontSize', 14, 'FontWeight', 'bold');
title('Temperature Profile', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
box on;

% Plot Upper Bound Conductivities
figure('Color', 'white', 'Units', 'inches', 'Position', [1, 1, 6, 11]);
subplot(2,1,1);
plot(depth, y_sigma_plus_1, '-', 'LineWidth', 2, 'DisplayName', '1 wt%');
hold on;
plot(depth, y_sigma_plus_01, '-', 'LineWidth', 2, 'DisplayName', '0.5 wt%');
plot(depth, y_sigma_plus_001, '-', 'LineWidth', 2, 'DisplayName', '0.1 wt%');
plot(depth, y_sigma_plus_0001, '-', 'LineWidth', 2, 'DisplayName', '1e-5 wt%');
hold off;
set(gca, 'YScale', 'log', 'FontSize', 12, 'LineWidth', 1.5);
ylim([1e-3, 10]);
xlim([200, 700]);
% xlabel('Depth (km)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Conductivity (S/m)', 'FontSize', 14, 'FontWeight', 'bold');
title('Yoshino model - HS upper bound', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
box on;
subplot(2,1,2);
plot(depth, k_sigma_plus_1, '-', 'LineWidth', 2, 'DisplayName', '1 wt%');
hold on;
plot(depth, k_sigma_plus_01, '-', 'LineWidth', 2, 'DisplayName', '0.1 wt%');
plot(depth, k_sigma_plus_001, '-', 'LineWidth', 2, 'DisplayName', '0.01 wt%');
plot(depth, k_sigma_plus_0001, '-', 'LineWidth', 2, 'DisplayName', '1e-5 wt%');
hold off;
set(gca, 'YScale', 'log', 'FontSize', 12, 'LineWidth', 1.5);
ylim([1e-3, 10]);
xlim([200, 700]);
xlabel('Depth (km)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Conductivity (S/m)', 'FontSize', 14, 'FontWeight', 'bold');
title('Karato model - HS upper bound', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
box on;

% Export the plot to a high-quality PNG file
exportgraphics(gcf, 'Yoshino_AREPS_2013_fig8_ub.png', 'Resolution', 300);

% Plot Lower Bound Conductivities
figure('Color', 'white', 'Units', 'inches', 'Position', [1, 1, 6, 11]);
subplot(2,1,1);
plot(depth, y_sigma_minus_1, '-', 'LineWidth', 2, 'DisplayName', '1 wt%');
hold on;
plot(depth, y_sigma_minus_01, '-', 'LineWidth', 2, 'DisplayName', '0.5 wt%');
plot(depth, y_sigma_minus_001, '-', 'LineWidth', 2, 'DisplayName', '0.1 wt%');
plot(depth, y_sigma_minus_0001, '-', 'LineWidth', 2, 'DisplayName', '1e-5 wt%');
hold off;
set(gca, 'YScale', 'log', 'FontSize', 12, 'LineWidth', 1.5);
ylim([1e-3, 10]);
xlim([200, 700]);
% xlabel('Depth (km)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Conductivity (S/m)', 'FontSize', 14, 'FontWeight', 'bold');
title('Lower Bound Conductivity', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
box on;
subplot(2,1,2);
plot(depth, k_sigma_minus_1, '-', 'LineWidth', 2, 'DisplayName', '1 wt%');
hold on;
plot(depth, k_sigma_minus_01, '-', 'LineWidth', 2, 'DisplayName', '0.1 wt%');
plot(depth, k_sigma_minus_001, '-', 'LineWidth', 2, 'DisplayName', '0.01 wt%');
plot(depth, k_sigma_minus_0001, '-', 'LineWidth', 2, 'DisplayName', '1e-5 wt%');
hold off;
set(gca, 'YScale', 'log', 'FontSize', 12, 'LineWidth', 1.5);
ylim([1e-3, 10]);
xlim([200, 700]);
xlabel('Depth (km)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Conductivity (S/m)', 'FontSize', 14, 'FontWeight', 'bold');
title('Lower Bound Conductivity', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
box on;
% Export the plot to a high-quality PNG file
exportgraphics(gcf, 'Yoshino_AREPS_2013_fig8_lb.png', 'Resolution', 300);





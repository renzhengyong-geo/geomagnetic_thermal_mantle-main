clear all;
% Define the input 2D array for continental geotherm
% Format: [Temperature (°C), Pressure (kbar)]
continental_geotherm = [
    642.52, 13.793;
    735.17, 20.094;
    827.8, 25.946;
    920.43, 31.798;
    1013.06, 37.651;
    1105.71, 43.951;
    1198.37, 50.252;
    1282.23, 56.55;
    1353.03, 65.085;
    1410.71, 74.961;
    1419.45, 74.067;
    1433.26, 85.275;
    1451.48, 96.933;
    1465.28, 108.142;
    1483.48, 119.351;
    1497.26, 130.111
];
% oceanic_geotherm=[
%    590.22	,10.214;
% 614.13	,11.217;
% 647.83	,11.785;
% 681.52	,13.012;
% 718.48	,13.582;
% 752.17	,14.04;
% 781.52	,15.265;
% 820.65	,16.825;
% 848.91	,17.39;
% 856.52	,16.734;
% 883.7	,17.299;
% 908.7	,18.631;
% 935.87	,19.525;
% 960.87	,19.869;
% 1001.09	,21.869;
% 1029.35	,21.994;
% 1060.87	,22.561;
% 1100	,24.67;
% 1128.26	,25.785;
% 1165.22	,26.135;
% 1216.3	,27.811;
% 1260.87	,29.704;
% 1296.74	,31.592;
% 1320.65	,34.132;
% 1343.48	,36.343;
% 1352.17	,38.106;
% 1364.13	,42.068;
% 1368.48	,45.807;
% 1373.91	,51.634;
% 1380.43	,54.495;
% 1388.04	,58.015;
% 1391.3	,59.666;
% 1392.39	,64.062;
% 1394.57	,62.415;
% 1397.83	,68.021;
% 1405.43	,70.552;
% 1408.7	,73.301;
% 1413.04	,78.249;
% 1421.74	,84.298;
% 1427.17	,82.103;
% 1429.35	,89.467;
% 1435.87	,92.547;
% 1439.13	,95.736;
% 1443.48	,104.31;
% 1447.83	,100.356;
% 1452.17	,106.732;
% 1458.7	,114.867;
% 1466.3	,121.575;
% 1469.57	,117.731;
% 1473.91	,124.876;
% 1482.61	,128.287;
% 1486.96	,132.136;
% 1490.22	,135.764
% ];

geotherm = continental_geotherm;
% Convert temperature from °C to Kelvin
n_points = length(geotherm(:, 1));

% constants for the following references 
% 1. Keppler, H., & Bolfan-Casanova, N. (2006). Thermodynamics of water solubility and partitioning. 
%   Reviews in Mineralogy & Geochemistry, 193–230. https://doi.org/10.2138/rmg.2006.62.9
% 2. Mierdel, K., Keppler, H., Smyth, J. R., & Langenhorst, F. (2007). Water solubility in aluminous orthopyroxene 
%   and the origin of earth s asthenosphere. Science, 315(5810), 364–368. https://doi.org/10.1126/science.1135422
% 3. Dong, J., Fischer, R. A., Stixrude, L. P., & Lithgow‐Bertelloni, C. R. (2021). 
%    Constraining the Volume of Earth s Early Oceans With a Temperature‐Dependent Mantle Water 
%    Storage Capacity Model. AGU Advances, 2(1). https://doi.org/10.1029/2020av000323
A_orthopyroxene = 0.01354;         % ppm/bar for Al-free enstatite or orthopyroxene
Delta_H_orthopyroxene = -4563;     % J/mol for Al-free enstatite/orthopyroxene
Delta_V_orthopyroxene = 12.1;      % cm^3/mol for Al-free enstatite/orthopyroxene
A_orthopyroxene_Al = 0.042;        % ppm/bar^0.5 for Al-bearing orthopyroxene
Delta_H_orthopyroxene_Al = -79685; % J/mol for Al-bearing orthopyroxene
Delta_V_orthopyroxene_Al = 11.3;   % cm^3/mol for Al-bearing orthopyroxene
A_olivine = 0.0066;                % ppm/bar for olivine
Delta_V_olivine = 10.6;            % cm^3/mol for olivine
A_clinopyroxene =7.144;            % ppm/bar^0.5 for clinopyroxene
Delta_V_clinopyroxene= 8.019;      % cm^3/mol for clinopyroxene
A_garnet= 0.679;                   % ppm/bar^0.5 for garnet
Delta_V_garnet = 5.71;             % cm^3/mol for garnet
R = 8.314;                         % J/(mol·K), gas constant

% Initialize arrays for the solubility contributions
orthopyroxene_al_free_water_solubility = zeros(n_points,1);
orthopyroxene_al_bearing_water_solubility = zeros(n_points,1);
orthopyroxene_water_solubility = zeros(n_points,1);
olivine_water_solubility = zeros(n_points,1);
orthopyroxene_olivine_water_solubility= zeros(n_points,1);
clinopyroxene_water_solubility= zeros(n_points,1);
garnet_water_solubility= zeros(n_points,1);

water_partition_olivine_vs_orthopyroxene=zeros(n_points,1);
water_partition_olivine_vs_clinopyroxene=zeros(n_points,1);
water_partition_olivine_vs_garnet=zeros(n_points,1);

% Calculate the water solubility for each point
for i = 1:n_points
    T_Celsius = geotherm(i, 1); % Temperature in °C
    P_GPa =  geotherm(i, 2)*0.1;% Pressure in GPa, 1kbar=0.1Gpa
    
    % Calculate water fugacity using water_fugacity_functions
    f_H2O = water_fugacity_pitzer_sterner_1994('PSfugacity', P_GPa, T_Celsius); % Ensure PSfugacity expects pressure in GPa and temperature in °C

    % Convert fugacity from GPa to bars if needed
    f_H2O = f_H2O * 1e4; % Assuming the fugacity function returns value in GPa， 1Gpa=10kbar=10000bar

    % Al-free contribution to water solubility, cc/mol*Gpa=1000J/mol
    orthopyroxene_c_water_Al_free = A_orthopyroxene * f_H2O * exp((-Delta_H_orthopyroxene) / (R * (T_Celsius + 273.15))) * ...
              exp((-Delta_V_orthopyroxene * P_GPa * 1e3) / (R * (T_Celsius + 273.15)));
    
    % Al-bearing contribution to water solubility
    orthopyroxene_c_water_Al = A_orthopyroxene_Al * sqrt(f_H2O) * exp((-Delta_H_orthopyroxene_Al) / (R * (T_Celsius + 273.15))) * ...
                 exp((-Delta_V_orthopyroxene_Al * P_GPa * 1e3) / (R * (T_Celsius + 273.15)));
    
    % Olivine water solubility, H/10^6Si/Mpa=6*1e-3 ppm/bar
    olivine_c_water = A_olivine* f_H2O * exp((-Delta_V_olivine* P_GPa * 1e3) / (R * (T_Celsius + 273.15)));
    
    % Clinopyroxene water solubility
    clinopyroxene_c_water = A_clinopyroxene* sqrt(f_H2O) * exp((-Delta_V_clinopyroxene* P_GPa * 1e3) / (R * (T_Celsius + 273.15)));

    % Garnet water solubility
    garnet_c_water = A_garnet* sqrt(f_H2O) * exp((-Delta_V_garnet* P_GPa * 1e3) / (R * (T_Celsius + 273.15)));

    % Store the individual contributions and total solubility
    orthopyroxene_al_free_water_solubility(i) = orthopyroxene_c_water_Al_free;
    orthopyroxene_al_bearing_water_solubility(i) = orthopyroxene_c_water_Al;
    orthopyroxene_water_solubility(i) = orthopyroxene_c_water_Al_free + orthopyroxene_c_water_Al;
    olivine_water_solubility(i) = olivine_c_water;
    orthopyroxene_olivine_water_solubility(i)=0.6*olivine_water_solubility(i)+orthopyroxene_water_solubility(i)*0.4;
    clinopyroxene_water_solubility(i) = clinopyroxene_c_water;
    garnet_water_solubility(i) = garnet_c_water;
    water_partition_olivine_vs_orthopyroxene(i)= olivine_water_solubility(i)/orthopyroxene_water_solubility(i);
    water_partition_olivine_vs_clinopyroxene(i) = olivine_water_solubility(i)/clinopyroxene_water_solubility(i);
    water_partition_olivine_vs_garnet(i)=olivine_water_solubility(i)/garnet_water_solubility(i);
end

% Create a high-quality figure for publication for water solubility
figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.35, 0.7]); % Adjusted figure size
hold on;

% Plot orthopyroxene with AL+ AL-free
plot(orthopyroxene_water_solubility, geotherm(:, 2), '-', 'LineWidth', 2.5, 'MarkerSize', 10,  'DisplayName', 'Al-saturated Orthopyroxene');
% Plot orthopyroxene Al-free
plot(orthopyroxene_al_free_water_solubility, geotherm(:, 2), '--b', 'LineWidth', 2.5, 'MarkerSize', 10, 'DisplayName', 'Al-free Orthopyroxene');
% Plot olivine 
plot(olivine_water_solubility, geotherm(:, 2), '-', 'LineWidth', 2.5, 'MarkerSize', 10, 'DisplayName', 'Olivine');
% Plot orthopyroxene+olivine
plot(orthopyroxene_olivine_water_solubility, geotherm(:, 2), '-r', 'LineWidth', 2.5, 'MarkerSize', 10,  'DisplayName', 'Al-Orthopyroxene*40%+Olivine*60%');
% Plot clinopyroxene
plot(clinopyroxene_water_solubility, geotherm(:, 2), '-', 'LineWidth', 2.5, 'MarkerSize', 10,  'DisplayName', 'Clinopyroxene');
% Plot garnet
plot(garnet_water_solubility, geotherm(:, 2), '-', 'LineWidth', 2.5, 'MarkerSize', 10,  'DisplayName', 'Garnet');

% Define the range for the low-velocity zone (between 40 kbar and 60 kbar)
y_min = 35; % lower bound for low-velocity zone in kbar
y_max = 60; % upper bound for low-velocity zone in kbar

% Create the shaded rectangle between 40 kbar and 60 kbar on the y-axis and 0 to 3000 ppm on the x-axis
fill([0, 3000, 3000, 0], [y_min, y_min, y_max, y_max], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Low-Velocity Zone');

% High-quality settings for the plot
xlabel('Water Solubility (ppm)', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');
ylabel('Pressure (kbar)', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');
grid off;
box on;
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontName', 'Arial'); % Set axis properties for better visualization
xlim([0, 3000]); % Set x-axis limits from 0 to 3000 ppm
ylim([0, 200]); % Set y-axis limits from 0 to 140 kbar
legend('FontSize', 14, 'Location', 'south', 'FontName', 'Arial');

% Reverse the y-axis direction
set(gca, 'YDir', 'reverse');

% Move the x-axis to the top and y-axis to the left
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'left');

% Improve grid appearance
%grid minor;
%set(gca, 'GridColor', [0.5, 0.5, 0.5], 'MinorGridColor', [0.8, 0.8, 0.8]); % Set grid line colors

% Save the figure in high-quality format (300 dpi)
saveas(gcf, 'water_solubility_land.png');
print(gcf, 'water_solubility_land', '-dpng', '-r300'); % Save at 300 dpi
hold off;


% Create a high-quality figure for publication for water partition
figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.35, 0.7]); % Adjusted figure size
hold on;

plot(water_partition_olivine_vs_orthopyroxene, geotherm(:, 2), '-', 'LineWidth', 2.5, 'MarkerSize', 10,  'DisplayName', 'Dolivine/orthopyroxene');
plot(water_partition_olivine_vs_clinopyroxene, geotherm(:, 2), '--b', 'LineWidth', 2.5, 'MarkerSize', 10, 'DisplayName', 'Dolivine/clinopyroxene');
plot(water_partition_olivine_vs_garnet, geotherm(:, 2), '-', 'LineWidth', 2.5, 'MarkerSize', 10, 'DisplayName', 'Dolivine/garnet/');

% Define the range for the low-velocity zone (between 40 kbar and 60 kbar)
y_min = 35; % lower bound for low-velocity zone in kbar
y_max = 60; % upper bound for low-velocity zone in kbar

% Create the shaded rectangle between 40 kbar and 60 kbar on the y-axis and 0 to 3000 ppm on the x-axis
fill([0, 3000, 3000, 0], [y_min, y_min, y_max, y_max], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Low-Velocity Zone');

% High-quality settings for the plot
xlabel('Water partition', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');
ylabel('Pressure (kbar)', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');
grid off;
box on;
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontName', 'Arial'); % Set axis properties for better visualization
% set(gca, 'XScale', 'log')
xlim([0, 1.5]); % Set x-axis limits from 0 to 
ylim([0, 200]); % Set y-axis limits from 0 to 200 kbar
legend('FontSize', 14, 'Location', 'south', 'FontName', 'Arial');

% Reverse the y-axis direction
set(gca, 'YDir', 'reverse');

% Move the x-axis to the top and y-axis to the left
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'left');

% Improve grid appearance
%grid minor;
%set(gca, 'GridColor', [0.5, 0.5, 0.5], 'MinorGridColor', [0.8, 0.8, 0.8]); % Set grid line colors

% Save the figure in high-quality format (300 dpi)
saveas(gcf, 'water_partition_land.png');
print(gcf, 'water_partition_land', '-dpng', '-r300'); % Save at 300 dpi
hold off;



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

% Convert temperature from °C to Kelvin
n_points = length(continental_geotherm(:, 1));

% Define constants for Al-free and Al-bearing water solubility
A = 0.01354; % ppm/bar for Al-free enstatite
A_Al = 0.042; % ppm/bar^0.5 for Al-bearing enstatite
Delta_H_1bar = -4563; % J/mol for Al-free enstatite
Delta_H_1bar_Al = -79685; % J/mol for Al-bearing enstatite
Delta_V_solid = 12.1; % cm^3/mol for Al-free enstatite
Delta_V_solid_Al = 11.3; % cm^3/mol for Al-bearing enstatite
R = 8.314; % J/(mol·K), gas constant

% Initialize arrays for the solubility contributions
total_water_solubility = zeros(n_points,1);
al_free_water_solubility = zeros(n_points,1);
al_bearing_water_solubility = zeros(n_points,1);
olivine_water_solubility = zeros(n_points,1);
upper_mantle_water_solubility= zeros(n_points,1);

% Calculate the water solubility for each point
for i = 1:n_points
    T_Celsius = continental_geotherm(i, 1); % Temperature in °C
    P_GPa =  continental_geotherm(i, 2)*0.1;% Pressure in GPa, 1kbar=0.1Gpa
    
    % Calculate water fugacity using water_fugacity_functions
    f_H2O = water_fugacity_pitzer_sterner_1994('PSfugacity', P_GPa, T_Celsius); % Ensure PSfugacity expects pressure in GPa and temperature in °C

    % Convert fugacity from GPa to bars if needed
    f_H2O = f_H2O * 1e4; % Assuming the fugacity function returns value in GPa， 1Gpa=10kbar=10000bar

    % Al-free contribution to water solubility, cc/mol*Gpa=1000J/mol
    c_water = A * f_H2O * exp((-Delta_H_1bar) / (R * (T_Celsius + 273.15))) * ...
              exp((-Delta_V_solid * P_GPa * 1e3) / (R * (T_Celsius + 273.15)));
    
    % Al-bearing contribution to water solubility
    c_water_Al = A_Al * sqrt(f_H2O) * exp((-Delta_H_1bar_Al) / (R * (T_Celsius + 273.15))) * ...
                 exp((-Delta_V_solid_Al * P_GPa * 1e3) / (R * (T_Celsius + 273.15)));
    
    % Olivine water solubility, H/10^6Si/Mpa=6*1e-3 ppm/bar
    c_water_ol = 1.1*(6*1e-3)* f_H2O * exp((-10.6 * P_GPa * 1e3) / (R * (T_Celsius + 273.15)));
    
    % Store the individual contributions and total solubility
    al_free_water_solubility(i) = c_water;
    al_bearing_water_solubility(i) = c_water_Al;
    total_water_solubility(i) = c_water + c_water_Al;
    olivine_water_solubility(i) = c_water_ol;
    upper_mantle_water_solubility(i)=0.6*olivine_water_solubility(i)+total_water_solubility(i)*0.4;
end

% Create a high-quality figure for publication
figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.35, 0.7]); % Adjusted figure size
hold on;

% Plot total water solubility (both contributions)
plot(total_water_solubility, continental_geotherm(:, 2), '-k', 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Al-bearing Orthopyroxene');

% Plot Al-free contribution
plot(al_free_water_solubility, continental_geotherm(:, 2), '--b', 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Al-free Orthopyroxene');

% Plot Al-bearing contribution
% plot(al_bearing_water_solubility, continental_geotherm(:, 2), '-^', 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Al-bearing Orthopyroxene');

% Plot olivine contribution
plot(olivine_water_solubility, continental_geotherm(:, 2), '--g', 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Olivine');

% Plot upper mantle contribution
plot(upper_mantle_water_solubility, continental_geotherm(:, 2), '-r', 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Al-Orthopyroxene*40%+Olivine*60%');


% Define the range for the low-velocity zone (between 40 kbar and 60 kbar)
y_min = 35; % lower bound for low-velocity zone in kbar
y_max = 60; % upper bound for low-velocity zone in kbar

% Create the shaded rectangle between 40 kbar and 60 kbar on the y-axis and 0 to 3000 ppm on the x-axis
fill([0, 3000, 3000, 0], [y_min, y_min, y_max, y_max], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Low-Velocity Zone');

% High-quality settings for the plot
xlabel('Water Solubility (ppm)', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');
ylabel('Pressure (kbar)', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');
grid on;
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontName', 'Arial'); % Set axis properties for better visualization
xlim([0, 3000]); % Set x-axis limits from 0 to 3000 ppm
ylim([0, 200]); % Set y-axis limits from 0 to 140 kbar
legend('FontSize', 14, 'Location', 'best', 'FontName', 'Arial');

% Reverse the y-axis direction
set(gca, 'YDir', 'reverse');

% Move the x-axis to the top and y-axis to the left
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'left');

% Improve grid appearance
grid minor;
set(gca, 'GridColor', [0.5, 0.5, 0.5], 'MinorGridColor', [0.8, 0.8, 0.8]); % Set grid line colors

% Save the figure in high-quality format (300 dpi)
saveas(gcf, 'water_solubility_contributions_vs_pressure_kbar.png');
print(gcf, 'water_solubility_contributions_vs_pressure_kbar', '-dpng', '-r300'); % Save at 300 dpi
hold off;


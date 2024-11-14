% Plot water fugacity vs. pressure for four different temperatures

% Define constants and coefficients
coeff = zeros(10, 6);
coeff(1, :) = [0, 0, 0.24657688e6, 0.51359951e2, 0, 0];
coeff(2, :) = [0, 0, 0.58638965e0, -0.28646939e-2, 0.31375577e-4, 0];
coeff(3, :) = [0, 0, -0.62783840e1, 0.14791599e-1, 0.35779579e-3, 0.15432925e-7];
coeff(4, :) = [0, 0, 0, -0.42719875e0, -0.16325155e-4, 0];
coeff(5, :) = [0, 0, 0.56654978e4, -0.16580167e2, 0.76560762e-1, 0];
coeff(6, :) = [0, 0, 0, 0.10917883e0, 0, 0];
coeff(7, :) = [0.38878656e13, -0.13494878e9, 0.30916564e6, 0.75591105e1, 0, 0];
coeff(8, :) = [0, 0, -0.65537898e5, 0.18810675e3, 0, 0];
coeff(9, :) = [-0.14182435e14, 0.18165390e9, -0.19769068e6, -0.23530318e2, 0, 0];
coeff(10, :) = [0, 0, 0.92093375e5, 0.12246777e3, 0, 0];

% Gas constant
R = 8314510; % Pa*cc/(K*mol)

% Temperatures in Kelvin to be considered
temperatures = [1273, 1473, 1673, 1873]; % Kelvin

% Pressures to evaluate (in GPa)
pressure_range = linspace(0.1, 15, 50); % Pressure from 0.1 GPa to 5 GPa

% Initialize figure
figure;
hold on;

% Colors for the plot lines
colors = ['b', 'r', 'g', 'm'];

% Loop through each temperature and calculate fugacity for each pressure
for t = 1:length(temperatures)
    temperature_K = temperatures(t);
    fugacity_values = zeros(size(pressure_range));
    
    for p = 1:length(pressure_range)
        pressure_GPa = pressure_range(p);
        fugacity_values(p) = PSfugacity(pressure_GPa, temperature_K - 273.15, coeff, R);
    end
    
    % Plot fugacity vs. pressure
    plot(pressure_range, fugacity_values, 'Color', colors(t), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('T = %d K', temperature_K));
end

% Set y-axis to logarithmic scale
set(gca, 'YScale', 'log');
% Set y-axis limits to [10^-5, 10^10]
ylim([1e-5, 1e10]);

% Labels and legend
xlabel('Pressure (GPa)', 'FontSize', 12);
ylabel('Fugacity (GPa)', 'FontSize', 12);
title('Water Fugacity vs Pressure for Different Temperatures', 'FontSize', 14);
legend('show', 'Location', 'NorthWest');
grid on;
box on; 
hold off;

% Save the figure
saveas(gcf, 'water_fugacity_vs_pressure.png'); % Saves the figure as a PNG file

% Function to calculate pressure using EOS
function pressure = PSeos(volume, temperature_C, targetP_GPa, coeff, R)
    % Convert target pressure from GPa to bars
    targetP = targetP_GPa * 1e4; % GPa to bars (1 GPa = 10,000 bars)
    % Convert temperature from °C to Kelvin
    temperature = temperature_C + 273.15;

    % Calculate density
    den = 1 ./ volume; % mol/cc
    
    % Calculate temperature-dependent coefficients
    c = zeros(1, 10);
    for i = 1:10
        c(i) = coeff(i, 1) * temperature^(-4) + coeff(i, 2) * temperature^(-2) + ...
               coeff(i, 3) * temperature^(-1) + coeff(i, 4) + ...
               coeff(i, 5) * temperature + coeff(i, 6) * temperature^2;
    end

    % Calculate pressure based on Pitzer and Sterner EOS
    numerator = (c(3) + 2 * c(4) * den + 3 * c(5) * den.^2 + 4 * c(6) * den.^3);
    denominator = (c(2) + c(3) * den + c(4) * den.^2 + c(5) * den.^3 + c(6) * den.^4).^2;
    pressure = (den + c(1) * den.^2 - den.^2 .* (numerator ./ denominator) + ...
                c(7) * den.^2 .* exp(-c(8) * den) + c(9) * den.^2 .* exp(-c(10) * den)) ...
                .* (R * temperature) / 1e5; % convert to bars

    % Return the difference between calculated pressure and target pressure
    pressure = pressure - targetP; % Difference used for root finding
end

% Function to calculate molar volume for a given pressure and temperature
function volume = PSvolume(pressure_GPa, temperature_C, coeff, R)
    % Use fminsearch to find the molar volume that gives the target pressure
    options = optimset('Display', 'off', 'TolX', 1e-8);
    % Starting guess for volume, optimized with bounds for realistic values
    initial_guess = 10; % Initial guess for molar volume in cc/mol
    % Objective function that finds the root of the pressure difference
    volume_obj = @(v) abs(PSeos(v, temperature_C, pressure_GPa, coeff, R));
    volume = fminsearch(volume_obj, initial_guess, options);
    
    % Ensure positive volume
    if volume <= 0
        error('Volume calculation failed: volume is non-positive.');
    end
    
    % Calculate pressure to check accuracy
    calculated_pressure_diff = PSeos(volume, temperature_C, pressure_GPa, coeff, R);
    tolerance = 1e-3; % Allowable tolerance for pressure difference (bars)
    if abs(calculated_pressure_diff) > tolerance
        error('Accurate molar volume not found: pressure difference exceeds tolerance.');
    end
end

% Function to calculate fugacity for a given pressure and temperature
function fugacity = PSfugacity(pressure_GPa, temperature_C, coeff, R)
    % Calculate the molar volume
    volume = PSvolume(pressure_GPa, temperature_C, coeff, R);
    
    % Convert temperature from °C to Kelvin
    temperature = temperature_C + 273.15;

    % Calculate density
    den = 1 / volume; % mol/cc

    % Calculate temperature-dependent coefficients
    c = zeros(1, 10);
    for i = 1:10
        c(i) = coeff(i, 1) * temperature^(-4) + coeff(i, 2) * temperature^(-2) + ...
               coeff(i, 3) * temperature^(-1) + coeff(i, 4) + ...
               coeff(i, 5) * temperature + coeff(i, 6) * temperature^2;
    end

    % Fugacity calculation
    fug = exp(log(den) + c(1) * den + ...
              (1 / (c(2) + c(3) * den + c(4) * den^2 + c(5) * den^3 + c(6) * den^4) - 1 / c(2)) - ...
              c(7) / c(8) * (exp(-c(8) * den) - 1) - ...
              c(9) / c(10) * (exp(-c(10) * den) - 1) + ...
              pressure_GPa * 1e4 * 1e5 / (den * R * temperature) + ...
              log(R * temperature) - 1) / 1e5; % convert to GPa

    fugacity = fug / 1e4; % Convert fugacity from bars to GPa (1 GPa = 10,000 bars)
end

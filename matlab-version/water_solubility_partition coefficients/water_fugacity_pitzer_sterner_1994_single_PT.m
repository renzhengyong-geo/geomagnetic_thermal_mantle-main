% Main Script to Use Pitzer and Sterner EOS
% Example usage
clear; clc;

% Coefficients for the EOS
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

pressure_GPa = 10; % Input Pressure in GPa
temperature_C = 10000; % Input Temperature in °C

try
    % Calculate molar volume and fugacity
    vol = PSvolume(pressure_GPa, temperature_C, coeff, R);
    fug = PSfugacity(pressure_GPa, temperature_C, coeff, R);
    
    fprintf('Molar Volume: %.6f cc/mol\n', vol);
    fprintf('Fugacity: %.6f GPa\n', fug);
catch ME
    % Display error message if accurate volume could not be found
    fprintf('Error: %s\n', ME.message);
end

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

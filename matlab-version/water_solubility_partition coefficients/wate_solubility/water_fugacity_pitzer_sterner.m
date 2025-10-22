function f_H2O = water_fugacity_pitzer_sterner(P, T)
% Calculate water fugacity using Pitzer & Sterner (1994) equation of state
% Input:
%   P - Pressure in GPa
%   T - Temperature in °C
% Output:
%   f_H2O - Water fugacity in GPa

    % Define coefficients for H2O from Pitzer & Sterner (1994) Table I
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

    % Gas constant: Pa*cc/(K*mol)
    R = 8314510;
    
    % Calculate molar volume first
    V = calculate_molar_volume(P, T, coeff, R);
    
    % Convert temperature to Kelvin
    T_K = T + 273.15;
    
    % Calculate density (mol/cc)
    rho = 1 / V;
    
    % Calculate temperature-dependent coefficients
    c = calculate_temperature_coefficients(T_K, coeff);
    
    % Convert pressure to Pa for calculation
    P_Pa = P * 1e9;  % GPa → Pa
    
    % Calculate compression factor z
    z = P_Pa / (rho * R * T_K);
    
    % Calculate residual Helmholtz energy (A_res/RT) from Pitzer & Sterner Eq.(1)
    A_res_RT = c(1) * rho + ...
               (1/(c(2) + c(3)*rho + c(4)*rho^2 + c(5)*rho^3 + c(6)*rho^4) - 1/c(2)) - ...
               c(7)/c(8) * (exp(-c(8)*rho) - 1) - ...
               c(9)/c(10) * (exp(-c(10)*rho) - 1);
    
    % Calculate fugacity coefficient
    ln_phi = A_res_RT + (z - 1) - log(z);
    phi = exp(ln_phi);
    
    % Fugacity in GPa
    f_H2O = phi * P;
end

function V = calculate_molar_volume(P, T, coeff, R)
% Calculate molar volume using fminsearch
    options = optimset('Display', 'off', 'TolX', 1e-10);
    
    if P > 1
        initial_guess = 5;
    else
        initial_guess = 20;
    end
    
    volume_obj = @(v) abs(calculate_pressure_difference(v, T, P, coeff, R));
    V = fminsearch(volume_obj, initial_guess, options);
    
    if V <= 0
        error('Volume calculation failed: non-positive volume found.');
    end
end

function pressure_diff = calculate_pressure_difference(V, T, P_target, coeff, R)
    T_K = T + 273.15;
    rho = 1 / V;
    
    c = calculate_temperature_coefficients(T_K, coeff);
    
    numerator = c(3) + 2*c(4)*rho + 3*c(5)*rho^2 + 4*c(6)*rho^3;
    denominator = (c(2) + c(3)*rho + c(4)*rho^2 + c(5)*rho^3 + c(6)*rho^4)^2;
    
    P_calc = (rho + c(1)*rho^2 - rho^2 * (numerator ./ denominator) + ...
              c(7)*rho^2 * exp(-c(8)*rho) + c(9)*rho^2 * exp(-c(10)*rho)) * ...
              (R * T_K) / 1e5;  % Convert to bar
    
    P_target_bar = P_target * 1e4;  % GPa → bar
    pressure_diff = P_calc - P_target_bar;
end

function c = calculate_temperature_coefficients(T, coeff)
    c = zeros(1, 10);
    for i = 1:10
        c(i) = coeff(i, 1) * T^(-4) + coeff(i, 2) * T^(-2) + ...
               coeff(i, 3) * T^(-1) + coeff(i, 4) + ...
               coeff(i, 5) * T + coeff(i, 6) * T^2;
    end
end
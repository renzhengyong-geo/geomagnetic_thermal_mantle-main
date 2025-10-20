function [C_OH, C_H2O_ppm] = olivine_water_solubility(P_GPa, T_C)
% Calculate water solubility in olivine using Kohlstedt et al. (1996) model

    % Kohlstedt et al. (1996) parameters
    A = 1.1;                    % H/10^6 Si per MPa at 1100°C
    n = 1;                      % Fugacity exponent
    delta_V = 10.6e-6;          % m³/mol
    R = 8.3145;                 % J/(mol·K)
    
    % Convert temperature to Kelvin
    T_K = T_C + 273.15;
    
    % Calculate water fugacity
    f_H2O_GPa = water_fugacity_PitzerSterner(P_GPa, T_C);
    
    % Convert to MPa for solubility calculation
    f_H2O_MPa = f_H2O_GPa * 1000;
    
    % Convert pressure to MPa
    P_MPa = P_GPa * 1000;
    
    % Calculate the exponential term - FIXED: ensure positive values
    exponent = -P_MPa * delta_V / (R * T_K);
    
    % Calculate hydroxyl solubility
    C_OH = A * (f_H2O_MPa^n) * exp(exponent);
    
    % Convert to wt ppm H2O
    C_H2O_ppm = C_OH * 0.06;  % Simplified conversion
end
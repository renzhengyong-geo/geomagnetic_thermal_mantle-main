function c_water_cpx_ppm = cpx_water_solubility(P, T)
% CPX_WATER_SOLUBILITY Calculate water solubility in chromian diopside
% based on the general solubility law from Bromiley et al. (2004)
%
% The function implements the general solubility equation:
% c_OH = A * fH2O^n * exp(-(P*ΔV)/(R*T))
%
% Equation Source:
%   Bromiley, G. D., Keppler, H., McCammon, C., Bromiley, F. A., & 
%   Jacobsen, S. D. (2004). Hydrogen solubility and speciation in natural,
%   gem-quality chromian diopside. American Mineralogist, 89(7), 941-949.
%   https://doi.org/10.2138/am-2004-0703
%   Equation 1
%
% Parameters for sample DI-1 (Chromian Diopside):
%   A = 2.15 ppm/bar^0.5, n = 0.5, ΔV = 7.43 cm³/mol
%
% Input:
%   P - Pressure in GPa
%   T - Temperature in °C
%
% Output:
%   c_water_cpx - H2O solubility in chromian diopside (ppm by weight)

    % Physical constants
    R = 8.3145;          % J/(mol·K) - Universal gas constant
    
    % ---------------------------------------------------------------------
    % Parameters for water solubility in chromian diopside from Bromiley et al. (2004)
    % These parameters are specific to natural Cr-diopside sample DI-1
    % ---------------------------------------------------------------------
    A = 2.15;            % Pre-exponential factor (ppm/bar^0.5)
    n = 0.5;             % Fugacity exponent (square root dependence)
    dV = 7.43;           % Volume change for incorporation (cm³/mol)
    
    % Convert temperature from Celsius to Kelvin for thermodynamic calculations
    T_K = T + 273.15;
    
    % Calculate water fugacity using the Pitzer & Sterner equation of state
    % Note: water_fugacity_pitzer_sterner returns fugacity in GPa
    fH2O_GPa = water_fugacity_pitzer_sterner(P, T);
    
    % Convert pressure and fugacity from GPa to bar
    % Critical unit conversion: 1 GPa = 10,000 bar
    P_bar = P * 1e4;           % Convert pressure from GPa to bar
    fH2O_bar = fH2O_GPa * 1e4; % Convert fugacity from GPa to bar
    
    % Convert volume change from cm³/mol to J/(mol·bar) for unit consistency
    % Conversion factor: 1 cm³·bar = 0.1 J
    % Therefore: ΔV (J/(mol·bar)) = ΔV (cm³/mol) × 0.1
    dV_bar = dV * 0.1;  % Convert cm³/mol to J/(mol·bar)
    
    % ---------------------------------------------------------------------
    % Calculate water solubility in chromian diopside
    % Using the general solubility law: c_OH = A * fH2O^n * exp(-(P*ΔV)/(R*T))
    %
    % Unit analysis:
    % - A: ppm/bar^0.5
    % - fH2O_bar: bar
    % - fH2O_bar^n: bar^0.5
    % - A * fH2O_bar^n: ppm/bar^0.5 * bar^0.5 = ppm
    % - Exponent: dimensionless [bar * J/(mol·bar) / (J/(mol·K) * K) = 1]
    % - Final result: ppm
    % ---------------------------------------------------------------------
    exponent = -(P_bar * dV_bar) / (R * T_K);
    c_water_cpx_ppm = A * (fH2O_bar)^n * exp(exponent);
end
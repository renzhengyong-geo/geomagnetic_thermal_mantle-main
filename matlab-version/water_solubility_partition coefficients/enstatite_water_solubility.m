function [c_water_total, c_water_Al_free, c_water_Al] = enstatite_water_solubility(P, T)
% CALCULATE_ENSTATITE_WATER_SOLUBILITY Calculate water solubility in enstatite
% according to the thermodynamic model of Mierdel et al. (2007)
%
% Reference:
%   Mierdel, K., Keppler, H., Smyth, J. R., & Langenhorst, F. (2007). 
%   Water solubility in aluminous orthopyroxene and the origin of Earth's 
%   asthenosphere. Science, 315(5810), 364-368. 
%   https://doi.org/10.1126/science.1135422
%
% Model Description:
%   This function calculates water solubility in orthopyroxene (enstatite) 
%   considering two independent incorporation mechanisms:
%   1. Water solubility in pure, Al-free enstatite (Eq. 1 in paper)
%   2. Additional water solubility coupled with Aluminum substitution (Eq. 2 in paper)
%   The total water solubility is the sum of these two components.
%
% Input:
%   P - Pressure in GPa
%   T - Temperature in °C
%
% Output:
%   c_water_total    - Total H2O solubility in aluminous enstatite (ppm by weight)
%   c_water_Al_free  - H2O solubility in pure, Al-free enstatite (ppm by weight)
%   c_water_Al       - Additional H2O solubility from Al-coupling (ppm by weight)
%
% Note: This function requires water_fugacity_pitzer_sterner(P, T) for 
%       water fugacity calculations.

    % Physical constants with consistent units
    R = 8.3145;          % J/(mol·K) = (kg·m²)/(s²·mol·K) - Universal gas constant
    
    % ---------------------------------------------------------------------
    % Parameters for water solubility in Al-free enstatite 
    % (Equation 1 from Mierdel et al., 2007)
    % Incorporation mechanism: 2H+ ↔ Mg2+ (forming OH pairs)
    % ---------------------------------------------------------------------
    A_Al_free = 0.01354; % Pre-exponential factor (ppm/bar)
    dH_Al_free = -4563;  % Enthalpy change for incorporation (J/mol)
    dV_Al_free = 12.1;   % Volume change for incorporation (cm³/mol)
    
    % ---------------------------------------------------------------------
    % Parameters for Al-coupled water solubility
    % (Equation 2 from Mierdel et al., 2007)
    % Incorporation mechanism: Al3+ + H+ ↔ Si4+ (forming isolated OH groups)
    % Note the square-root dependence on water fugacity
    % ---------------------------------------------------------------------
    A_Al = 0.042;        % Pre-exponential factor (ppm/bar^0.5)
    dH_Al = -79685;      % Enthalpy change for Al-coupled incorporation (J/mol)
    dV_Al = 11.3;        % Volume change for Al-coupled incorporation (cm³/mol)
    
    % Convert temperature from Celsius to Kelvin for thermodynamic calculations
    T_K = T + 273.15;
    
    % Calculate water fugacity using the Pitzer & Sterner (1994) equation of state
    % Reference: Pitzer, K. S., & Sterner, S. M. (1994). J. Chem. Phys.
    fH2O_GPa = water_fugacity_pitzer_sterner(P, T);
    
    % Convert pressure and fugacity from GPa to bar
    % Note: The original equations in Mierdel et al. (2007) use bar units
    % 1 GPa = 10,000 bar
    P_bar = P * 1e4;           % Convert pressure from GPa to bar
    fH2O_bar = fH2O_GPa * 1e4; % Convert fugacity from GPa to bar
    
    % Convert volume changes from cm³/mol to J/(mol·bar) for unit consistency
    % 1 cm³ = 1e-6 m³, 1 bar = 1e5 Pa, 1 J = 1 Pa·m³
    % Therefore: 1 cm³·bar = 1e-6 m³ · 1e5 Pa = 0.1 J
    % So: ΔV (J/(mol·bar)) = ΔV (cm³/mol) × 0.1
    dV_Al_free_bar = dV_Al_free * 0.1;  % Convert cm³/mol to J/(mol·bar)
    dV_Al_bar = dV_Al * 0.1;            % Convert cm³/mol to J/(mol·bar)
    
    % ---------------------------------------------------------------------
    % Calculate water solubility in Al-free enstatite (Eq. 1)
    % c_water = A * fH2O * exp(-(ΔH + P*ΔV)/(R*T))
    % Units: ppm/bar * bar * exp(-(J/mol + bar * J/(mol·bar))/(J/(mol·K)*K))
    % All units cancel properly in the exponent
    % ---------------------------------------------------------------------
    exponent_Al_free = -(dH_Al_free + P_bar * dV_Al_free_bar) / (R * T_K);
    c_water_Al_free = A_Al_free * fH2O_bar * exp(exponent_Al_free);
    
    % ---------------------------------------------------------------------
    % Calculate additional water solubility coupled to Aluminum (Eq. 2)
    % c_water_Al = A_Al * sqrt(fH2O) * exp(-(ΔH_Al + P*ΔV_Al)/(R*T))
    % Units: ppm/bar^0.5 * bar^0.5 * exp(-(J/mol + bar * J/(mol·bar))/(J/(mol·K)*K))
    % All units cancel properly in the exponent
    % ---------------------------------------------------------------------
    exponent_Al = -(dH_Al + P_bar * dV_Al_bar) / (R * T_K);
    c_water_Al = A_Al * sqrt(fH2O_bar) * exp(exponent_Al);
    
    % ---------------------------------------------------------------------
    % Calculate total water solubility in aluminous enstatite
    % Total solubility is the sum of Al-free and Al-coupled components
    % This represents water solubility in orthopyroxene saturated with
    % an aluminous phase (spinel or garnet) under upper mantle conditions
    % ---------------------------------------------------------------------
    c_water_total = c_water_Al_free + c_water_Al;
end
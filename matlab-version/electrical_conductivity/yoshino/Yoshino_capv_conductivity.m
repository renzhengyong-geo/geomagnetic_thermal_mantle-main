function sigma = Yoshino_capv_conductivity(T, Cw, P)
    % ------------------------------------------------------------------------
    % FUNCTION: Yoshino_ca_pv_conductivity
    % ------------------------------------------------------------------------
    % DESCRIPTION:
    % Compute the electrical conductivity of Ca-perovskite (CaSiO3) based on
    % experimental data from Xu et al. (1998).
    %
    % This function models electrical conductivity dominated by ionic conduction
    % in Ca-perovskite, with negligible contributions from hopping and
    % proton conduction mechanisms.
    %
    % REFERENCES:
    % Xu Y, McCammon C, Poe B T. The effect of alumina on the electrical 
    % conductivity of silicate perovskite. Science, 1998, 282(5390): 922-924.
    %
    % Yoshino, T., & Katsura, T. (2013). Electrical conductivity of mantle minerals: 
    % Role of water in conductivity anomalies. Annual Review of Earth and Planetary Sciences, 
    % 41, 605–628.
    %
    % INPUTS (SCALARS ONLY):
    %   T   - Temperature in Kelvin
    %   Cw  - Water content in wt% (not used in current model)
    %   P   - Pressure in Pascals (currently not used - placeholder)
    %
    % OUTPUT:
    %   sigma - Electrical conductivity in S/m
    % ------------------------------------------------------------------------

    % Constants - all in J units
    k = 1.380649e-23;                    % Boltzmann constant, J/K
    eV_to_J = 1.60218e-19;               % Conversion factor from eV to J
    
    % Ca-perovskite parameters from Xu et al. (1998)
    % Converted to Joules for consistency
    sigma_0i  = 631;                     % S/m (Pre-exponential factor for ionic conduction)
    delta_Hi  = 1.03 * eV_to_J;          % J (Activation enthalpy for ionic conduction)
    
    % Hopping conduction negligible for Ca-perovskite
    sigma_0h  = 0;                       % S/m
    delta_Hh  = 0;                       % J
    
    % Proton conduction not significant in dry Ca-perovskite
    sigma_0p  = 0;                       % S/m  
    delta_H0_p = 0;                      % J
    alpha_p    = 0;                      % J

    % Input validation (scalar inputs)
    if T <= 0
        error('Temperature (T) must be greater than zero Kelvin.');
    end
    if Cw < 0 || Cw > 1
        error('Water content (Cw) must be in the range of 0 to 1 (wt%).');
    end
    if P < 0
        error('Pressure (P) must be non-negative.');
    end

    % Note: Current model based on dry Ca-perovskite measurements
    % Water content parameter is included for consistency but not used
    if Cw > 0
        warning('Water content effects on Ca-perovskite conductivity are not well constrained in this model');
    end

    % Calculate conduction mechanisms using J units
    % Dominated by ionic conduction in Ca-perovskite
    sigma_ionic = sigma_0i * exp(-delta_Hi / (k * T));
    
    % Hopping conduction (negligible)
    sigma_hop = sigma_0h * exp(-delta_Hh / (k * T));
    
    % Proton conduction (negligible in dry conditions)
    sigma_proton = sigma_0p * Cw * exp(-(delta_H0_p - alpha_p * Cw^(1/3)) / (k * T));
    
    % Total conductivity (dominated by ionic conduction)
    sigma = sigma_ionic + sigma_hop + sigma_proton;
end
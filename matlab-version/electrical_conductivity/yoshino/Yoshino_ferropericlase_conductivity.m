function sigma = Yoshino_ferropericlase_conductivity(T, Cw, P, Xfe)
    % ------------------------------------------------------------------------
    % FUNCTION: Yoshino_ferropericlase_conductivity
    % ------------------------------------------------------------------------
    % DESCRIPTION:
    % Compute the electrical conductivity of ferropericlase using Yoshino's 
    % experimental data with pressure and iron content dependence.
    %
    % This function models electrical conductivity dominated by small polaron
    % (hopping) conduction in ferropericlase, with explicit pressure dependence
    % and iron concentration effects.
    %
    % REFERENCES:
    % Yoshino T, Ito E, Katsura T, et al. Effect of iron content on electrical 
    % conductivity of ferropericlase with implications for the spin transition 
    % pressure. Journal of Geophysical Research: Solid Earth, 2011, 116(B4).
    %
    % Yoshino, T., & Katsura, T. (2013). Electrical conductivity of mantle minerals: 
    % Role of water in conductivity anomalies. Annual Review of Earth and Planetary Sciences, 
    % 41, 605–628.
    %
    % INPUTS (SCALARS ONLY):
    %   T   - Temperature in Kelvin
    %   Cw  - Water content in wt% (not used in current model)
    %   P   - Pressure in Pascals
    %   Xfe - Iron content as mole fraction (0 to 1)
    %         Typical mantle value: Xfe = 0.1 (10% iron content)
    %         From Yoshino & Katsura (2013), Page 6: "for mantle peridotite (X_Fe = 0.1)"
    %
    % OUTPUT:
    %   sigma - Electrical conductivity in S/m
    % ------------------------------------------------------------------------

    % Constants - all in J units
    k = 1.380649e-23;                    % Boltzmann constant, J/K
    eV_to_J = 1.60218e-19;               % Conversion factor from eV to J
    
    % Ferropericlase parameters from Yoshino et al. (2011)
    % Converted to Joules for consistency
    sigma_0i  = 0;                       % S/m (Ionic conduction negligible)
    delta_Hi  = 0;                       % J
    
    % Hopping conduction parameters
    sigma_0h  = 19;                      % S/m (Pre-exponential factor)
    delta_Hh0 = 0.63 * eV_to_J;          % J (Base activation enthalpy)
    alpha_h   = 0.66 * eV_to_J;          % J (Iron concentration dependence)
    
    % Pressure dependence parameters
    % Note: 0.0104 factor converts (GPa * cm³/mol) to eV
    % V0 is in cm³/mol, pressure is converted to GPa internally
    V0 = -0.45;                          % cm³/mol (Activation volume)
    beta_h = -0.61 * eV_to_J;            % J (Pressure-iron coupling parameter)
    
    % Proton conduction negligible for ferropericlase
    sigma_0p  = 0;                       % S/m
    delta_H0_p = 0;                      % J
    alpha_p   = 0;                       % J

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
    if Xfe < 0 || Xfe > 1
        error('Iron content (Xfe) must be in the range of 0 to 1.');
    end

    % Note: Current model based on dry ferropericlase measurements
    % Water content parameter is included for consistency but not used
    if Cw > 0
        warning('Water content effects on ferropericlase conductivity are not well constrained in this model');
    end

    % Convert pressure from Pa to GPa for activation volume calculation
    % 1 GPa = 1e9 Pa
    P_GPa = P * 1e-9;
    
    % Calculate conduction mechanisms using J units
    % Ionic conduction (negligible)
    sigma_ionic = sigma_0i * exp(-delta_Hi / (k * T));
    
    % Hopping conduction with pressure and iron content dependence
    % The 0.0104 factor converts (GPa * cm³/mol) to eV
    % We need to work in eV for the pressure term to maintain consistency
    % with the original formulation, then convert final activation energy to J
    
    % Calculate activation enthalpy in eV (for pressure term compatibility)
    delta_Hh_eV = (delta_Hh0 / eV_to_J) - alpha_h * Xfe^(1/3) / eV_to_J;
    
    % Add pressure contribution in eV (0.0104 converts GPa*cm³/mol to eV)
    pressure_term_eV = 0.0104 * P_GPa * (V0 - (beta_h / eV_to_J) * Xfe^(1/3));
    
    % Total activation enthalpy in eV, then convert to J
    delta_Hh_total_eV = delta_Hh_eV + pressure_term_eV;
    delta_Hh_total_J = delta_Hh_total_eV * eV_to_J;
    
    sigma_hop = sigma_0h * Xfe * exp(-delta_Hh_total_J / (k * T));
    
    % Proton conduction (negligible in dry conditions)
    sigma_proton = sigma_0p * Cw * exp(-(delta_H0_p - alpha_p * Cw^(1/3)) / (k * T));
    
    % Total conductivity (dominated by hopping conduction)
    sigma = sigma_ionic + sigma_hop + sigma_proton;
end
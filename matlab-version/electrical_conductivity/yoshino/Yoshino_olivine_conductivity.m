function sigma = Yoshino_olivine_conductivity(T, Cw, P, Xfe)
    % YOSHINO_OLIVINE_CONDUCTIVITY Calculate electrical conductivity of olivine
    %
    % References:
    %   Yoshino, T., & Katsura, T. (2013). Electrical conductivity of mantle minerals: 
    %   Role of water in conductivity anomalies. Annual Review of Earth and Planetary Sciences, 
    %   41, 605-628. https://doi.org/10.1146/annurev-earth-050212-124022
    %
    %   Yoshino, T., Matsuzaki, T., Shatskiy, A., & Katsura, T. (2009). The effect of 
    %   water on the electrical conductivity of olivine aggregates and its implications 
    %   for the electrical structure of the upper mantle. Earth and Planetary Science Letters, 
    %   288(1-2), 291-300. https://doi.org/10.1016/j.epsl.2009.09.032
    %
    % Inputs:
    %   T    - Temperature in Kelvin (K)
    %   Cw   - Water content in weight percent (wt%)
    %   P    - Pressure in Pa (not used in current model)
    %   Xfe  - Iron content (dimensionless, 0 to 1)
    %
    % Output:
    %   sigma - Electrical conductivity in S/m
    %
    % Model parameters for olivine:
    %   sigma_0i = 10^4.73 S/m, delta_Hi = 2.31 eV  (ionic conduction)
    %   sigma_0h = 10^2.98 S/m, delta_Hh = 1.71 eV  (small polaron conduction)  
    %   sigma_0p = 10^1.90 S/m, delta_H0_p = 0.92 eV (proton conduction)
    %   alpha_p = 0.16 eV/wt%^(1/3) (proton conduction parameter)

    % Constants
    k = 1.380649e-23;                    % Boltzmann constant [J/K]
    J2eV = 1.0 / 1.60217733e-19;         % Conversion factor: 1 eV = 1.60217733e-19 J
    k_eV = k * J2eV;                     % Boltzmann constant [eV/K]
    
    % Model parameters for olivine (Yoshino et al., 2009, 2013)
    sigma_0i = 10^4.73;                  % Pre-exponential factor for ionic conduction [S/m]
    delta_Hi = 2.31;                     % Activation enthalpy for ionic conduction [eV]
    
    sigma_0h = 10^2.98;                  % Pre-exponential factor for small polaron conduction [S/m]
    delta_Hh = 1.71;                     % Activation enthalpy for small polaron conduction [eV]
    
    sigma_0p = 10^1.90;                  % Pre-exponential factor for proton conduction [S/m]
    delta_H0_p = 0.92;                   % Base activation enthalpy for proton conduction [eV]
    alpha_p = 0.16;                      % Proton concentration dependence parameter [eV/wt%^(1/3)]
    
    % Input validation
    if T <= 0
        error('Temperature (T) must be greater than zero Kelvin.');
    end
    
    if Cw < 0 || Cw > 100
        error('Water content (Cw) must be in the range of 0 to 100 wt%%.');
    end
    
    if Xfe < 0 || Xfe > 1
        error('Iron content (Xfe) must be in the range of 0 to 1.');
    end
    
    if P < 0
        error('Pressure (P) must be non-negative.');
    end
    
    % Calculate individual conduction mechanisms
    ionic_conduction = sigma_0i * exp(-delta_Hi / (k_eV * T));
    polaron_conduction = sigma_0h * Xfe * exp(-delta_Hh / (k_eV * T));
    proton_conduction = sigma_0p * Cw * exp(-(delta_H0_p - alpha_p * Cw^(1/3)) / (k_eV * T));
    
    % Total conductivity (sum of all mechanisms)
    sigma = ionic_conduction + polaron_conduction + proton_conduction;
end
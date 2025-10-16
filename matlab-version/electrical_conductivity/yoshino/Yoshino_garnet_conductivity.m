function sigma = Yoshino_garnet_conductivity(T, Cw, P)
    % YOSHINO_GARNET_CONDUCTIVITY Calculate electrical conductivity of majorite garnet
    %
    % Reference:
    %   Yoshino T, Nishi M, Matsuzaki T, et al. Electrical conductivity of majorite 
    %   garnet and its implications for electrical structure in the mantle transition 
    %   zone. Physics of the Earth and Planetary Interiors, 2008, 170(3-4): 193-200.
    %
    % Inputs:
    %   T    - Temperature in Kelvin (K)
    %   Cw   - Water content in weight percent (wt%)
    %   P    - Pressure in Pa (not used in current model)
    %
    % Output:
    %   sigma - Electrical conductivity in S/m
    %
    % Model parameters for garnet:
    %   sigma_h = 1070 S/m, Hh = 1.59 eV (small polaron conduction)
    %   Note: Ionic conduction and proton conduction are negligible

    % Constants
    k = 1.380649e-23;                    % Boltzmann constant [J/K]
    J2eV = 1.0 / 1.60217733e-19;         % Conversion factor: 1 eV = 1.60217733e-19 J
    k_eV = k * J2eV;                     % Boltzmann constant [eV/K]
    
    % Model parameters for garnet (Yoshino et al., 2008)
    sigma_i = 0;                         % Ionic conduction (negligible for garnet)
    Hi = 0;                              % Activation enthalpy for ionic conduction [eV]
    
    sigma_h = 1070;                      % Pre-exponential factor for small polaron conduction [S/m]
    Hh = 1.59;                           % Activation enthalpy for small polaron conduction [eV]
    
    sigma_p = 0;                         % Proton conduction (negligible for garnet)
    Hp = 0;                              % Activation enthalpy for proton conduction [eV]
    alpha = 0;                           % Proton concentration dependence parameter [eV/wt%^(1/3)]
    
    % Input validation
    if T <= 0
        error('Temperature (T) must be greater than zero Kelvin.');
    end
    
    if Cw < 0 || Cw > 100
        error('Water content (Cw) must be in the range of 0 to 100 wt%%.');
    end
    
    if P < 0
        error('Pressure (P) must be non-negative.');
    end
    
    % Calculate individual conduction mechanisms
    ionic_conduction = sigma_i * exp(-Hi / (k_eV * T));
    polaron_conduction = sigma_h * exp(-Hh / (k_eV * T));
    proton_conduction = sigma_p * Cw * exp(-(Hp - alpha * Cw^(1/3)) / (k_eV * T));
    
    % Total conductivity (sum of all mechanisms)
    sigma = ionic_conduction + polaron_conduction + proton_conduction;
end
function sigma = Yoshino_clinopyroxene_conductivity(T, Cw, P)
    % YOSHINO_CLINOPYROXENE_CONDUCTIVITY Calculate electrical conductivity of clinopyroxene
    %
    % Reference:
    %   Xu Y, Shankland T J. Electrical conductivity of orthopyroxene and its 
    %   high pressure phases. Geophysical Research Letters, 1999, 26(17): 2645-2648.
    %
    % Inputs:
    %   T    - Temperature in Kelvin (K)
    %   Cw   - Water content in weight percent (wt%)
    %   P    - Pressure in Pa (not used in current model)
    %
    % Output:
    %   sigma - Electrical conductivity in S/m
    %
    % Model parameters for clinopyroxene:
    %   sigma_i = 1778 S/m, Hi = 1.87 eV (ionic conduction)
    %   Note: Proton conduction and small polaron conduction are negligible

    % Constants
    k = 1.380649e-23;                    % Boltzmann constant [J/K]
    J2eV = 1.0 / 1.60217733e-19;         % Conversion factor: 1 eV = 1.60217733e-19 J
    k_eV = k * J2eV;                     % Boltzmann constant [eV/K]
    
    % Model parameters for clinopyroxene (Xu & Shankland, 1999)
    sigma_i = 1778;                      % Pre-exponential factor for ionic conduction [S/m]
    Hi = 1.87;                           % Activation enthalpy for ionic conduction [eV]
    
    sigma_h = 0;                         % Small polaron conduction (negligible for clinopyroxene)
    Hh = 0;                              % Activation enthalpy for small polaron conduction [eV]
    
    sigma_p = 0;                         % Proton conduction (negligible for clinopyroxene)
    Hp = 0;                              % Activation enthalpy for proton conduction [eV]
    r = 0;                               % Water content exponent
    
    % Input validation
    if T <= 0
        error('Temperature (T) must be greater than zero Kelvin.');
    end
    
    if Cw < 0 || Cw > 200
        error('Water content (Cw) must be in the range of 0 to 200 wt%%.');
    end
    
    if P < 0
        error('Pressure (P) must be non-negative.');
    end
    
    % Calculate individual conduction mechanisms
    ionic_conduction = sigma_i * exp(-Hi / (k_eV * T));
    polaron_conduction = sigma_h * exp(-Hh / (k_eV * T));
    proton_conduction = sigma_p * Cw^r * exp(-Hp / (k_eV * T));
    
    % Total conductivity (sum of all mechanisms)
    sigma = ionic_conduction + polaron_conduction + proton_conduction;
end
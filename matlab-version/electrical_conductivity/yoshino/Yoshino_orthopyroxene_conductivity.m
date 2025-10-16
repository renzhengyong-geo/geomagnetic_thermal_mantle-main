function sigma = Yoshino_orthopyroxene_conductivity(T, Cw, P)
    % YOSHINO_ORTHOPYROXENE_CONDUCTIVITY Calculate electrical conductivity of orthopyroxene
    %
    % Reference:
    %   Zhang B, Yoshino T, Wu X, et al. Electrical conductivity of enstatite 
    %   as a function of water content: Implications for the electrical structure 
    %   in the upper mantle. Earth and Planetary Science Letters, 2012, 357: 11-20.
    %
    % Inputs:
    %   T    - Temperature in Kelvin 
    %   Cw   - Water content in weight percent (0-100 wt%)
    %   P    - Pressure in Pa 
    %
    % Output:
    %   sigma - Electrical conductivity in S/m
    %
    % Model parameters for orthopyroxene:
    %   sigma_i = 10^3.99 S/m, Hi = 1.88 eV  (ionic conduction)
    %   sigma_p = 10^2.58 S/m, Hp = 0.84 eV  (proton conduction)
    %   alpha = 0.08 eV/wt%^(1/3) (proton conduction parameter)

    % Constants
    k = 1.380649e-23;                    % Boltzmann constant [J/K]
    J2eV = 1.0 / 1.60217733e-19;         % Conversion factor: 1 eV = 1.60217733e-19 J
    k_eV = k * J2eV;                     % Boltzmann constant [eV/K]
    
    % Model parameters for orthopyroxene (Zhang et al., 2012)
    sigma_i = 10^3.99;                   % Pre-exponential factor for ionic conduction [S/m]
    Hi = 1.88;                           % Activation enthalpy for ionic conduction [eV]
    
    sigma_h = 0;                         % Small polaron conduction (negligible for orthopyroxene)
    Hh = 0;                              % Activation enthalpy for small polaron conduction [eV]
    
    sigma_p = 10^2.58;                   % Pre-exponential factor for proton conduction [S/m]
    Hp = 0.84;                           % Base activation enthalpy for proton conduction [eV]
    alpha = 0.08;                        % Proton concentration dependence parameter [eV/wt%^(1/3)]
    
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
    
    % Proton conduction with water content dependence
    proton_conduction = sigma_p * Cw * exp(-(Hp - alpha * Cw^(1/3)) / (k_eV * T));
    
    % Total conductivity (sum of all mechanisms)
    sigma = ionic_conduction + polaron_conduction + proton_conduction;
end

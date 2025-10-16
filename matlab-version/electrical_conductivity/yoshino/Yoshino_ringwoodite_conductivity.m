function sigma = Yoshino_ringwoodite_conductivity(T, Cw, P, Xfe)
    % ------------------------------------------------------------------------
    % FUNCTION: Yoshino_ringwoodite_conductivity
    % ------------------------------------------------------------------------
    % DESCRIPTION:
    % Compute the electrical conductivity of Ringwoodite using Yoshino's group 
    % high-temperature and high-pressure experiments.
    %
    % FORMULA:
    % The conductivity is computed using the formula:
    % sigma = sigma_0i * exp(-delta_Hi / (k * T)) + ...
    %         sigma_0h * Xfe * exp(-delta_Hh / (k * T)) + ...
    %         sigma_0p * Cw  * exp(-(delta_H0_p - alpha_p * Cw^(1/3)) / (k * T));
    %
    % REFERENCES:
    % Yoshino, T., & Katsura, T. (2013). Electrical conductivity of mantle minerals: 
    % Role of water in conductivity anomalies. Annual Review of Earth and Planetary Sciences, 
    % 41, 605–628. https://doi.org/10.1146/annurev-earth-050212-124022
    %
    % Yoshino T, Manthilake G, Matsuzaki T, Katsura T. 2008a. Dry mantle transition zone 
    % inferred from electrical conductivity of wadsleyite and ringwoodite. Nature 451:326–29
    %
    % INPUTS (SCALARS ONLY):
    %   T   - Temperature in Kelvin
    %   Cw  - Water content in wt% (typical range: 0.001 to 0.1 wt%)
    %   P   - Pressure in Pascals (currently not used - placeholder)
    %   Xfe - Iron content as mole fraction (0 to 1)
    %
    % OUTPUT:
    %   sigma - Electrical conductivity in S/m
    % ------------------------------------------------------------------------

    % Constants
    k = 1.380649e-23;           % Boltzmann constant, J/K
    eV_to_J = 1.60218e-19;      % Conversion factor
    
    % Ringwoodite parameters from Yoshino et al. 2008a (Table 1, Page 12)
    sigma_0i  = 0;                      % S/m (Ionic conduction negligible)
    delta_Hi  = 0;                      % J
    
    sigma_0h  = 10^(2.92);              % S/m (Hopping conduction pre-exponential)
    delta_Hh  = 1.36 * eV_to_J;         % J (Activation enthalpy for hopping)
    
    sigma_0p  = 10^(1.44);              % S/m (Proton conduction pre-exponential)
    delta_H0_p = 1.12 * eV_to_J;        % J (Base activation enthalpy for protons)
    alpha_p    = 0.67 * eV_to_J;        % J (Water concentration dependence)

    % Input validation (scalar inputs)
    if T <= 0
        error('Temperature (T) must be greater than zero Kelvin.');
    end
    if Cw < 0 || Cw > 100
        error('Water content (Cw) must be in the range of 0 to 100 (wt%).');
    end
    if Xfe < 0 || Xfe > 1
        error('Iron content (Xfe) must be in the range of 0 to 1.');
    end

    % Calculate conduction mechanisms
    sigma_ionic = sigma_0i * exp(-delta_Hi / (k * T));
    sigma_hop = sigma_0h * Xfe * exp(-delta_Hh / (k * T));
    
    % CORRECTED: Proper sign for proton conduction
    sigma_proton = sigma_0p * Cw * exp(-(delta_H0_p - alpha_p * Cw^(1/3)) / (k * T));
    
    % Total conductivity
    sigma = sigma_ionic + sigma_hop + sigma_proton;
    
    % Placeholder for pressure correction
    sigma = sigma * P^(0); % Pressure effect not considered, placeholder for consistency
end
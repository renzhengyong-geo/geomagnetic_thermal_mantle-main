% ------------------------------------------------------------------------
% FUNCTION: Yoshino_ringwoodite_conductivity
% ------------------------------------------------------------------------
% DESCRIPTION:
% Compute the conductivity of Ringwoodite using Yoshino's group high-temperature
% and high-pressure experiments.
%
% FORMULA:
% The conductivity is computed using the formula:
% sigma = sigma_0i * exp(-delta_Hi / (k * T)) + ...
%         sigma_0h * Xfe * exp(-delta_Hh / (k * T)) + ...
%         sigma_0p * Cw  * exp((-delta_H0_p + alpha_p * Cw^(1/3)) / (k * T));
%
% where:
%   sigma       : Electrical conductivity
%   sigma_0i    : Pre-exponential factor for ionic conduction 
%   delta_Hi    : Activation enthalpy for ionic conduction (in joules)
%   k           : Boltzmann constant (1.380649e-23 J/K)
%   T           : Absolute temperature (in Kelvin)
%   sigma_0h    : Pre-exponential factor for hopping of electron holes
%   Xfe         : Iron mole fraction or concentration (unitless)
%   delta_Hh    : Standard activation enthalpy for hopping conduction (in joules)
%   sigma_0p    : Pre-exponential factor for proton (H) conduction
%   Cw          : Water content or concentration (unitless)
%   delta_H0_p  : Standard activation enthalpy for proton conduction (in joules)
%   alpha_p     : Parameter related to concentration of water or proton-conducting species
%
% REFERENCES:
% Yoshino, T., & Katsura, T. (2013). Electrical conductivity of mantle minerals: 
% Role of water in conductivity anomalies. Annual Review of Earth and Planetary Sciences, 
% 41, 605–628. https://doi.org/10.1146/annurev-earth-050212-124022
% This study suggest that XFE=1.0 as its effect is considered in the
% sigma_0H!
%
% Yoshino T, Manthilake G, Matsuzaki T, Katsura T. 2008a.Dry mantle transition zone 
% inferred from electrical conductivity of wadsleyite and ringwoodite. Nature 451:326–29
%
% INPUTS:
%   T   - Temperature in Kelvin
%   Cw  - Water content (wt%, 0 to 100)
%   P   - Pressure (in Pascals)
%   Xfe - Iron content (dimensionless, 0 to 1)

% OUTPUT:
%   sigma - Electrical conductivity (S/m)
% ------------------------------------------------------------------------

function sigma = Yoshino_ringwoodite_conductivity(T, Cw, P, Xfe)
    % Constants
    k = 1.380649e-23;           % Boltzmann constant, J/K
    sigma_0i  = 0;              % S/m (Pre-exponential factor for ionic conduction)
    delta_Hi  = 0;              % J (Converted from eV to J)
    sigma_0h  = 10^(2.92);      % S/m (Pre-exponential factor for hopping conduction)
    delta_Hh  = 1.36 * 1.60218e-19;   % J (Converted from eV to J)
    sigma_0p  = 10^(1.44);            % S/m (Pre-exponential factor for proton conduction)
    delta_H0_p = 1.12 * 1.60218e-19;  % J (Converted from eV to J)
    alpha_p    = 0.67 * 1.60218e-19;  % J (Converted from eV to J)

    % Input validation
    if any(T <= 0)
        error('Temperature (T) must be greater than zero Kelvin.');
    end
    if any(Cw < 0) || any(Cw > 100)
        error('Water content (Cw) must be in the range of 0 to 100 (wt%).');
    end
    if any(Xfe < 0) || any(Xfe > 1)
        error('Iron content (Xfe) must be in the range of 0 to 1.');
    end

    % Sum the conductivities
    sigma = sigma_0i * exp(-delta_Hi / (k * T)) + ...
            sigma_0h * Xfe * exp(-delta_Hh / (k * T)) + ...
            sigma_0p * Cw  * exp((-delta_H0_p + alpha_p * Cw^(1.0/3.0)) / (k * T));
    
    % Placeholder for pressure correction
    sigma = sigma * P^(0); % Pressure effect not considered, placeholder for consistency
end

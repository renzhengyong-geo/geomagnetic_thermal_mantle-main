% ------------------------------------------------------------------------
% FUNCTION: Yoshino_Olivine_conductivity
% ------------------------------------------------------------------------
% DESCRIPTION:
% Compute the conductivity of olivine using Yoshino's group high-temperature
% and high-pressure experiments.
%
% FORMULA:
% The conductivity is computed using the formula:
% sigma = sigma_0i * exp(-delta_Hi / (k * T)) + ...
%         sigma_0h * Xfe * exp(-delta_Hh / (k * T)) + ...
%         sigma_0p * Cw  * exp((-delta_H0_p + alpha_p * Cw^(1/3)) / (k * T));
%
% REFERENCES:
% Yoshino, T., & Katsura, T. (2013). Electrical conductivity of mantle minerals: 
% Role of water in conductivity anomalies. Annual Review of Earth and Planetary Sciences, 
% 41, 605–628. https://doi.org/10.1146/annurev-earth-050212-124022
%
% Yoshino, T., Matsuzaki, T., Shatskiy, A., & Katsura, T. (2009). The effect of 
% water on the electrical conductivity of olivine aggregates and its implications 
% for the electrical structure of the upper mantle. Earth and Planetary Science Letters, 
% 288(1–2), 291–300. https://doi.org/10.1016/j.epsl.2009.09.032
%
% INPUTS:
%   T   - Temperature in Kelvin
%   Xfe - Fe content (dimensionless, range: 0 to 1)
%   Cw  - Water content in wt% (range: 0 to 100)
%   P   - Pressure (arbitrary, kept for consistent structure, in Pascal)
%
% OUTPUT:
%   sigma - Electrical conductivity (S/m)
% ------------------------------------------------------------------------

function sigma = Yoshino_olivine_conductivity(T, Cw, P, Xfe)
    % Constants
    k = 1.380649e-23;             % Boltzmann constant, J/K
    sigma_0i  = 10^(4.73);        % S/m (Pre-exponential factor for ionic conduction)
    delta_Hi  = 2.31 * 1.60218e-19; % J (Converted from eV to J)
    sigma_0h  = 10^(2.98);        % S/m (Pre-exponential factor for hopping conduction)
    delta_Hh  = 1.71 * 1.60218e-19; % J (Converted from eV to J)
    sigma_0p  = 10^(1.90);        % S/m (Pre-exponential factor for proton conduction)
    delta_H0_p = 0.92 * 1.60218e-19; % J (Converted from eV to J)
    alpha_p   = 0.16 * 1.60218e-19; % J (Converted from eV to J)

    % Input validation
    if any(T <= 0)
        error('Temperature must be greater than zero Kelvin.');
    end

    if any(Cw < 0) || any(Cw > 100)
        error('Water content (Cw) must be in the range of 0 to 100 wt%.');
    end

    if any(Xfe < 0) || any(Xfe > 1)
        error('Iron content (Xfe) must be in the range of 0 to 1.');
    end

    % Conductivity calculation
    sigma = sigma_0i * exp(-delta_Hi / (k * T)) + ...
            sigma_0h * Xfe * exp(-delta_Hh / (k * T)) + ...
            sigma_0p * Cw  * exp((-delta_H0_p + alpha_p * Cw^(1.0/3.0)) / (k * T));

    % Pressure correction placeholder (P^0 is 1)
    sigma = sigma * P^(0); % Pressure effect not considered, placeholder for consistency
end

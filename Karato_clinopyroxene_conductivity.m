% ------------------------------------------------------------------------
% FUNCTION: Karato_orthopyroxene_conductivity
% ------------------------------------------------------------------------
% DESCRIPTION:
% Compute the conductivity of Orthopyroxene using Karato's group high temperature
% and pressure experiments.
%
% FORMULA:
% The conductivity is computed using the formula:
%   sigma = A1*Cw^{r1}*exp(-(E1+PV1)/(RT)) + A2*Cw^{r2}*exp(-(E2+PV2)/(RT))
% where:
%   - A1, A2 are the pre-exponential terms which already contain the influence
%     of the Mg# and oxygen fugacity.
%   - Cw is the water content.
%   - r1, r2 are constants that depend on the mechanism of electrical conduction.
%   - E1, E2 are the activation energies.
%   - V1, V2 are the activation volumes.
%   - Subscript 1 stands for polaron conduction.
%   - Subscript 2 stands for proton conduction (e.g., H).
%
% REFERENCES:
% Values of A, r, E, V are from:
%   Karato, S. I. (2011). Water distribution across the mantle transition zone
%   and its implications for global material circulation. Earth and Planetary
%   Science Letters, 301(3–4), 413–423. https://doi.org/10.1016/j.epsl.2010.11.038
%
% UNITS:
%   - temperature: Kelvin (K); 0 K is equivalent to -273.15°C.
%   - water_content: dimensionless, typically in the range 0-1.
%   - pressure: Pascal (Pa); 1 Pa = 1 J/m³ = 1e-6 J/cc.
%
% INPUTS:
%   temperature   - Temperature in Kelvin.
%   water_content - Water content (dimensionless).
%   pressure      - Pressure in Pascal.
%
% OUTPUT:
%   sigma - Electrical conductivity (S/m).
% ------------------------------------------------------------------------

function sigma = Karato_clinopyroxene_conductivity(temperature, water_content, pressure)
    % Constants
    R = 8.314;       % Gas constant, J/(K⋅mol)
    A1 = 10^(2.7);   % Pre-exponential factor for polaron conduction, S/m
    A2 = 10^(2.6);   % Pre-exponential factor for proton conduction, S/m
    r1 = 0;          % Empirical constant for polaron conduction
    r2 = 0.62;       % Empirical constant for proton conduction
    E1 = 147000;     % Activation energy for polaron conduction, J/mol
    E2 = 82000;      % Activation energy for proton conduction, J/mol
    V1 = 0.0;        % Activation volume for polaron conduction, cc/mol
    V2 = 0.0;        % Activation volume for proton conduction, cc/mol

    % Check if temperature is greater than zero
    if any(temperature <= 0)
        error('Temperature must be greater than zero Kelvin.');
    end

    % Ensure pressure is numeric
    if iscell(pressure)
        pressure = cell2mat(pressure); % Convert cell array to numeric if necessary
    end

    % Check if pressure is non-negative
    if any(pressure < 0)
        error('Pressure must be non-negative.');
    end

    % Convert pressure from Pa to J/cc
    pressure = 1e3 * pressure;  % 1 Pa = 1 J/m³ = 1e-6 J/cc

    % Ensure water_content is numeric
    if iscell(water_content)
        water_content = cell2mat(water_content); % Convert cell array to numeric if necessary
    end

    % Check if water_content is within the range of 0 to 1
    if any(water_content < 0) || any(water_content > 1)
        error('Water content must be in the range of 0 to 1.');
    end

    % Calculate conductivities
    sigma_polaron = A1 * water_content.^r1 .* exp(-(E1 + pressure * V1) / (R * temperature));
    sigma_proton =  A2 * water_content.^r2 .* exp(-(E2 + pressure * V2) / (R * temperature));

    % Sum the conductivities
    sigma = sigma_polaron + sigma_proton;
end
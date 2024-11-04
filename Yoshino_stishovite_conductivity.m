%Yoshino T, Shimojuku A, Li D. Electrical conductivity of stishovite as a function of water content[J]. 
%Physics of the Earth and Planetary Interiors, 2014, 227: 48-54.

function sigma = Yoshino_stishovite_conductivity(temperature, water_content,pressure)
    % Constants
    % This function calculates the electrical conductivity based on the
    % provided temperature (T), water content (Cw), and material constants.
    % Constants
    k = 1.380649*10^(-23);
    J2eV = 1.0/(1.60217733*10^(-19)); % 1eV = 1.60217733*10^(-19) J
    k = k*J2eV;
    
    sigma_i = 137;  % S/m
    Hi = 1.03;  % activation enthalpy, J/mol
    sigma_h = 0;  % S/m
    Hh = 0;  % S/m
    sigma_p = 1042;
    Hp = 1.07;  % S/m
    alpha = 0.19;
    
    T =  temperature;
    Cw = water_content;
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
    pressure = 1e-6 * pressure;  % 1 Pa = 1 J/m³ = 1e-6 J/cc

    % Ensure water_content is numeric
    if iscell(water_content)
        water_content = cell2mat(water_content); % Convert cell array to numeric if necessary
    end

    % Check if water_content is within the range of 0 to 1
    if any(water_content < 0) || any(water_content > 1)
        error('Water content must be in the range of 0 to 1.');
    end

    % Sum the conductivities
    sigma = sigma_i*exp(-Hi/(k*T))+...
        sigma_h*exp(-Hh/(k*T))+...
        sigma_p*Cw*exp(-(Hp-alpha*Cw^(1.0/3.0))/(k*T));
end

%Katsura T, Yokoshi S, Kawabe K, et al. Pressure dependence of electrical conductivity of (Mg, Fe) SiO 3 ilmenite[J]. 
%Physics and Chemistry of Minerals, 2007, 34: 249-255.

function sigma = Yoshino_akimotoite_conductivity(temperature, water_content,pressure)
    % Constants
    % This function calculates the electrical conductivity based on the
    % provided temperature (T), water content (Cw), and material constants.
    % Constants
    k = 1.380649*10^(-23);
    J2eV = 1.0/(1.60217733*10^(-19)); % 1eV = 1.60217733*10^(-19) J
    k = k*J2eV;
    sigma_i = 15  - 1;  % S/m
    Va = -1.5;
    Hi = (0.82 - 0.01) + pressure * Va * 0.0104;
    sigma_h = 0;  % S/m
    Hh = 0;  % S/m
    sigma_p = 0;
    Hp = 0;  % S/m
    r = 0;
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
        sigma_p*Cw^r*exp(-Hp/(k*T));
end

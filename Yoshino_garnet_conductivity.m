%Yoshino T, Nishi M, Matsuzaki T, et al. Electrical conductivity of majorite garnet and its implications for electrical structure in the mantle transition zone[J]. 
%Physics of the Earth and Planetary Interiors, 2008, 170(3-4): 193-200.
function sigma = Yoshino_garnet_conductivity(temperature, water_content,pressure)
    % Constants
    % This function calculates the electrical conductivity based on the
    % provided temperature (T), water content (Cw), and material constants.
    % Constants
    k = 1.380649*10^(-23);
    J2eV = 1.0/(1.60217733*10^(-19)); % 1eV = 1.60217733*10^(-19) J
    k = k*J2eV;
    T =  temperature;
    Cw = water_content;
    sigma_i = 0;
    Hi = 0;
    sigma_h = 1070;
    Hh = 1.59 ;
    %sigma_h = 522;% - 10;
    %Hh = 1.47 ;%+ 0.03;
    sigma_p = 0;
    Hp = 0;
    alpha = 0;     
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

    % Sum the conductivities
    sigma = sigma_i*exp(-Hi/(k*T))+...
        sigma_h*exp(-Hh/(k*T))+...
        sigma_p*Cw*exp(-(Hp-alpha*Cw^(1.0/3.0))/(k*T));
end

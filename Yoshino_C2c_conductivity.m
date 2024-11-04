%Zhang B, Yoshino T, Wu X, et al. Electrical conductivity of enstatite as a function of water content: Implications for the electrical structure in the upper mantle[J]. 
%Earth and Planetary Science Letters, 2012, 357: 11-20.


function sigma = Yoshino_C2c_conductivity(temperature, water_content,pressure)
    % Constants
    % This function calculates the electrical conductivity based on the
    % provided temperature (T), water content (Cw), and material constants.
    % Constants
    k = 1.380649*10^(-23);
    J2eV = 1.0/(1.60217733*10^(-19)); % 1eV = 1.60217733*10^(-19) J
    k = k*J2eV;
    
    sigma_i = 10^3.9;
    Hi = 1.88;
    sigma_h = 0;
    Hh = 0 ;
    sigma_p = 10^2.58 ;
    Hp = 0.84 ;
    alpha = 0.08 ;
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

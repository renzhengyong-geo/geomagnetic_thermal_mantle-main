%Yoshino T, Ito E, Katsura T, et al. Effect of iron content on electrical conductivity of ferropericlase with implications for the spin transition pressure[J]. 
%Journal of Geophysical Research: Solid Earth, 2011, 116(B4).

function sigma = Yoshino_ferropericlase_conductivity_Xfe(temperature, water_content,pressure,Xfe)
    % Constants
    % This function calculates the electrical conductivity based on the
    % provided temperature (T), water content (Cw), and material constants.
    % Constants
    k = 1.380649*10^(-23);
    J2eV = 1.0/(1.60217733*10^(-19)); % 1eV = 1.60217733*10^(-19) J
    k = k*J2eV;
    sigma_i = 0;  % S/m
    Hi = 0;
    
    sigma_h = 19 ;%240;%19 ;  % S/m
    Hh = 0.63 ;%1.08;%  ;  % S/m
    alpha_h = 0.66 ;
    belta_h = -0.61;
    V0 = -0.45 ;
    
    sigma_p = 0;
    Hp = 0;  % S/m
    alpha = 0;
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
        sigma_h*Xfe*exp(-((Hh-alpha_h*Xfe^(1.0/3.0) + 0.0104 * pressure * (V0 - belta_h * Xfe^(1.0/3.0))))/(k*T))+...
        sigma_p*Cw*exp(-(Hp-alpha*Cw^(1.0/3.0))/(k*T));
end

function fO2 = compute_fO2_Oneill(temperature_C)
    % COMPUTE_FO2_ONEILL computes oxygen fugacity (fO2) using the O'Neill (1986) method.
    %  1. O'Neill, H. St. C. (1986). Mo-MoO₂ (MOM) oxygen buffer and the free energy 
    %    of formation of MoO₂. American Mineralogist, 71, 1007–1010. 
    
    % Inputs:
    %   temperature_C - Temperature in Celsius (°C)
    %
    % Outputs:
    %   fO2 - Oxygen fugacity (fO2) in bar

    % Constants
    R = 8.314; % Universal gas constant, J/(mol·K)

    % Convert temperature from Celsius to Kelvin
    temperature_K = temperature_C + 273.15;

    % Chemical potential of oxygen for the O'Neill method
    mu_O2 = -603268 + 337.460 .* temperature_K - 20.6892 .* temperature_K .* log(temperature_K);

    % Calculate fO2
    fO2 = exp(mu_O2 ./ (R .* temperature_K));
end
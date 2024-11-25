function output = water_partition_wad_ring_gt(Temperature, Pressure)
    % WATER_PARTITION_WAD_RING_GT Calculates the partitioning of water between 
    % minerals wadsleyite, ringwoodite, and garnet in the mantle transition zone.
    %
    % INPUTS:
    %   Temperature - Temperature in degrees Celsius (°C)
    %   Pressure    - Pressure in kilobars (kbar)
    %
    % OUTPUT:
    %   output - Array containing water partition ratios:
    %            [wad_water_solubility, ring_water_solubility, Wadsleyite/Ringwoodite, Wadsleyite/Garnet]
    %
    % REFERENCES:
    %   1. Dong, J., Fischer, R. A., Stixrude, L. P., & Lithgow‐Bertelloni, C. R. (2021).
    %      Constraining the Volume of Earth’s Early Oceans With a Temperature‐Dependent Mantle Water
    %      Storage Capacity Model. AGU Advances, 2(1), e2020AV000323.
    %      https://doi.org/10.1029/2020AV000323
    %
    %   2. Inoue, T., Katsuda, M., Yurimoto, H., & Irifune, T. (2006). 
    %      H2O partitioning between wadsleyite, ringwoodite, and garnet.

    % === Convert temperature and pressure to appropriate units ===
    T_Kelvin = Temperature + 273.15; % Convert °C to Kelvin
    P_GPa = Pressure * 0.1;          % Convert kbar to GPa (1 kbar = 0.1 GPa)
    
    % === Constants for water solubility calculations ===
    % Regression coefficients for wadsleyite (wad) and ringwoodite (ring)
    % Solubility model: ln(c_H2O) = a + b / T
    wad_a = -7.6356;                 % Coefficient 'a' for wadsleyite
    wad_b = 13739.5371;              % Coefficient 'b' for wadsleyite
    ring_a = -6.8856;                % Coefficient 'a' for ringwoodite
    ring_b = 12206.2676;             % Coefficient 'b' for ringwoodite

    % === Water solubility calculations ===
    wad_water_solubility = exp(wad_a + wad_b / T_Kelvin); % Wadsleyite solubility
    ring_water_solubility = exp(ring_a + ring_b / T_Kelvin); % Ringwoodite solubility

    % === Calculate water partition ratios ===
    water_partition_wad_vs_ring = wad_water_solubility / ring_water_solubility; % Wad/Ring ratio
    water_partition_wad_vs_gt = 10; % Wad/Garnet ratio (constant, based on reference)

    % === Output array ===
    output = [
        wad_water_solubility, ...          % Solubility of wadsleyite
        ring_water_solubility, ...         % Solubility of ringwoodite
        water_partition_wad_vs_ring, ...   % Wad/Ring ratio
        water_partition_wad_vs_gt          % Wad/Garnet ratio
    ];
end

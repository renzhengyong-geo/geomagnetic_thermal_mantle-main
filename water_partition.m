function output = water_partition(Temperature, Pressure)
    % water_partition calculates the partitioning of water between minerals.
    %
    % INPUTS:
    %   Temperature - Temperature in degrees Celsius (°C)
    %   Pressure - Pressure in kilobars (kbar)
    %
    % OUTPUT:
    %   output - Array containing water partition ratios:
    %            [Olivine/Orthopyroxene, Olivine/Clinopyroxene, Olivine/Garnet, Wadsleyite/ringwoodite, Wadsleyite/Garnet]

    % === References ===
    % Constants for the following references:
    % 1. Keppler, H., & Bolfan-Casanova, N. (2006). Thermodynamics of water solubility and partitioning.
    %    Reviews in Mineralogy & Geochemistry, 62, 193–230.
    %    https://doi.org/10.2138/rmg.2006.62.9
    %
    % 2. Mierdel, K., Keppler, H., Smyth, J. R., & Langenhorst, F. (2007). Water solubility in aluminous orthopyroxene
    %    and the origin of Earth's asthenosphere. Science, 315(5810), 364–368.
    %    https://doi.org/10.1126/science.1135422
    %
    % 3. Dong, J., Fischer, R. A., Stixrude, L. P., & Lithgow‐Bertelloni, C. R. (2021).
    %    Constraining the Volume of Earth’s Early Oceans With a Temperature‐Dependent Mantle Water
    %    Storage Capacity Model. AGU Advances, 2(1), e2020AV000323.
    %    https://doi.org/10.1029/2020AV000323

    % 4. Toru Inoue; Masanori Katsuda; Hisayoshi Yurimoto; Tetsuo Irifune. (2006). 
    %    H2O partitioning between wadsleyite and garnet and between ringwoodite and garnet.

    % === Mineral constants ===
    % Orthopyroxene (Al-free)
    A_orthopyroxene = 0.01354;         % ppm/bar
    Delta_H_orthopyroxene = -4563;     % J/mol
    Delta_V_orthopyroxene = 12.1;      % cm^3/mol
    
    % Orthopyroxene (Al-bearing)
    A_orthopyroxene_Al = 0.042;        % ppm/bar^0.5
    Delta_H_orthopyroxene_Al = -79685; % J/mol
    Delta_V_orthopyroxene_Al = 11.3;   % cm^3/mol
    
    % Olivine
    A_olivine = 0.0066;                % ppm/bar
    Delta_V_olivine = 10.6;            % cm^3/mol
    
    % Clinopyroxene
    A_clinopyroxene = 7.144;           % ppm/bar^0.5
    Delta_V_clinopyroxene = 8.019;     % cm^3/mol
    
    % Garnet
    A_garnet = 0.679;                  % ppm/bar^0.5
    Delta_V_garnet = 5.71;             % cm^3/mol
    
    % Gas constant
    R = 8.314; % J/(mol·K)
    
    % === Temperature and pressure adjustments ===
    T_Kelvin = Temperature + 273.15;   % Convert °C to Kelvin
    P_GPa = Pressure * 0.1;            % Convert kbar to GPa (1 kbar = 0.1 GPa)
    
    % === Calculate water fugacity ===
    f_H2O = water_fugacity_pitzer_sterner_1994('PSfugacity', P_GPa, Temperature);
    f_H2O = f_H2O * 1e4;  % Convert fugacity from GPa to bars (1 GPa = 10,000 bars)
    
    % === Water solubility calculations ===
    % Al-free orthopyroxene
    orthopyroxene_c_water_Al_free = A_orthopyroxene * f_H2O * ...
        exp((-Delta_H_orthopyroxene) / (R * T_Kelvin)) * ...
        exp((-Delta_V_orthopyroxene * P_GPa * 1e3) / (R * T_Kelvin));
    
    % Al-bearing orthopyroxene
    orthopyroxene_c_water_Al = A_orthopyroxene_Al * sqrt(f_H2O) * ...
        exp((-Delta_H_orthopyroxene_Al) / (R * T_Kelvin)) * ...
        exp((-Delta_V_orthopyroxene_Al * P_GPa * 1e3) / (R * T_Kelvin));
    
    % Olivine
    olivine_c_water = A_olivine * f_H2O * ...
        exp((-Delta_V_olivine * P_GPa * 1e3) / (R * T_Kelvin));
    
    % Clinopyroxene
    clinopyroxene_c_water = A_clinopyroxene * sqrt(f_H2O) * ...
        exp((-Delta_V_clinopyroxene * P_GPa * 1e3) / (R * T_Kelvin));
    
    % Garnet
    garnet_c_water = A_garnet * sqrt(f_H2O) * ...
        exp((-Delta_V_garnet * P_GPa * 1e3) / (R * T_Kelvin));
    
    % % Regression coefficients are shown with standard errors (SE)
    % % Wad, a=-7.6356 (0.8952), b=13739.5371 (1602.2032)
    % % ringwoodite  a=-6.8856 (1.3651) b=12206.2676 (2400.7947)
    % % ln_c_H2O = a + b / T; % ln(c_H2O)
    % % c_H2O = exp(ln_c_H2O); % c_H2O
    % Wad = exp(a+b/T_Kelvin)

    % === Calculate total solubility and partition ratios ===
    orthopyroxene_water_solubility = orthopyroxene_c_water_Al_free + orthopyroxene_c_water_Al;
    olivine_water_solubility = olivine_c_water;
    clinopyroxene_water_solubility = clinopyroxene_c_water;
    garnet_water_solubility = garnet_c_water;
    
    % Water partition ratios
    water_partition_olivine_vs_orthopyroxene = olivine_water_solubility / orthopyroxene_water_solubility;
    water_partition_olivine_vs_clinopyroxene = olivine_water_solubility / clinopyroxene_water_solubility;
    water_partition_olivine_vs_garnet = olivine_water_solubility / garnet_water_solubility;
    
    % === Output ===
    output = [water_partition_olivine_vs_orthopyroxene, ...
              water_partition_olivine_vs_clinopyroxene, ...
              water_partition_olivine_vs_garnet];
end



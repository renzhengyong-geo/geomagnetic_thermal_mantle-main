function sigma = carbonate_melt_conductivity(T)
    % Function to calculate electrical conductivity for carbonate melts 
    % using the Arrhenius equation.
    % 
    % Citation:
    % Fabrice Gaillard; Mohammed Malki; Giada Iacono-Marziano; Michel Pichavant; 
    % Bruno Scaillet. (2008). Carbonatite Melts and Electrical Conductivity in 
    % the Asthenosphere. Science, 322(5906), 1360–1363. 
    % https://doi.org/10.1126/science.1163033
    %
    % Constants
    R = 8.314;  % Universal gas constant in J/mol·K
    
    % Activation energies (E_a) and pre-exponential factors (sigma_0) for each composition
    % The table below provides the values used in the Arrhenius equation for each composition:
    %
    % Composition                           Ea (J/mol)      sigma_0 (S/m)
    % ----------------------------------------------------------------------
    % (LiNaK)_2(CO_3)_3                    32,500          6,590
    % (NaK)_2(CO_3)_2                      31,427          4,177
    % (NaKCa_0.5)_2(CO_3)_3                32,500          4,144
    % (NaKCa)(CO_3)_2                      30,307          2,504
    % (KCa_0.5)_2(CO_3)_2                  34,489          3,149
    % Mantle Carbonatites                  31,900          3,440
    %
    % Ea is the activation energy in J/mol, and sigma_0 is the pre-exponential 
    % factor in S/m for each composition.
    %
    % Ensure T is a column vector to match dimensions for element-wise operations
    T = T(:);  % Convert T to a column vector
    
    % Activation energies (Ea) and pre-exponential factors (sigma_0)
    Ea = [32500, 31427, 32500, 30307, 34489, 31900];  % J/mol
    sigma_0 = [6590, 4177, 4144, 2504, 3149, 3440];     % S/m
    
    sigma = zeros(6,1);
    % Calculate conductivity for each composition using the Arrhenius equation
    % Ea and sigma_0 are now both 6x1 vectors, so they can broadcast across T
    for i=1:6
        sigma(i) = sigma_0(i) * exp(-Ea(i) / (R * T));  % Element-wise multiplication
    end

end
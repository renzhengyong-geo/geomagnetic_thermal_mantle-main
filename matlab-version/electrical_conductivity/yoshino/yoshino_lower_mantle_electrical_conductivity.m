function [sigma_upper, sigma_lower] = yoshino_lower_mantle_electrical_conductivity(T, P, pv_f, fp_f, capv_f, C_water, water_partition_coefficients)
    % ------------------------------------------------------------------------
    % FUNCTION: yoshino_lower_mantle_electrical_conductivity
    % ------------------------------------------------------------------------
    % DESCRIPTION:
    % Compute electrical conductivity of LOWER MANTLE using Yoshino's 
    % experimental data for silicate perovskite, ferropericlase, and 
    % calcium perovskite.
    %
    % The lower mantle (660-2900 km depth) is dominated by silicate perovskite
    % (bridgmanite), ferropericlase, and calcium perovskite according to 
    % pyrolitic mantle composition.
    %
    % INPUTS:
    %   T - Temperature in Kelvin
    %   P - Pressure in Pa  
    %   pv_f - Volume fraction of silicate perovskite (bridgmanite) (0-1)
    %   fp_f - Volume fraction of ferropericlase (0-1)
    %   capv_f - Volume fraction of calcium perovskite (0-1)
    %   C_water - Total water content in wt% (typical: 0.001-0.1 wt%)
    %   water_partition_coefficients - Partition coefficients [silicate_perovskite, ferropericlase, calcium_perovskite]
    %
    % REFERENCES:
    % Yoshino, T., & Katsura, T. (2013). Electrical conductivity of mantle minerals: 
    % Role of water in conductivity anomalies. Annual Review of Earth and Planetary Sciences, 
    % 41, 605–628.
    % ------------------------------------------------------------------------

    % Input validation
    if T <= 0
        error('Temperature T must be positive in K, but is %.6f', T);
    end
    if P < 0
        error('Pressure P must be non-negative in Pa, but is %.6f', P);
    end
    
    % Validate mineral fractions
    mineral_fractions = [pv_f, fp_f, capv_f];
    if any(mineral_fractions < 0)
        error('Mineral fractions must be positive: pv_f=%.6f, fp_f=%.6f, capv_f=%.6f', pv_f, fp_f, capv_f);
    end
    if abs(sum(mineral_fractions) - 1.0) > 1e-10
        error('Mineral fractions must sum to 1.0, but sum to %.10f', sum(mineral_fractions));
    end
    
    % Validate water content (typical mantle range)
    if C_water < 0
        error('C_water must be non-negative, but is %.6f', C_water);
    end
    if C_water > 1.0
        warning('Water content %.4f wt%% may be unrealistically high for mantle conditions', C_water);
    end
    
    % Validate partition coefficients
    if length(water_partition_coefficients) ~= 3
        error('water_partition_coefficients must have exactly 3 elements');
    end
    if abs(sum(water_partition_coefficients) - 1.0) > 1e-10
        error('water_partition_coefficients must sum to 1.0, but sum to %.10f', sum(water_partition_coefficients));
    end

    % Standard mantle composition: Xfe = Fe/(Fe+Mg) = 0.1
    % From Yoshino & Katsura (2013), Page 6: "for mantle peridotite (X_Fe = 0.1)"
    Xfe_in_ferropericlase = 0.1; 
    
    % Calculate individual mineral conductivities with partitioned water
    water_in_perovskite = C_water * water_partition_coefficients(1);
    water_in_ferropericlase = C_water * water_partition_coefficients(2);
    water_in_capv = C_water * water_partition_coefficients(3);    

    % Use Yoshino's conductivity models for lower mantle minerals
    sigma_pv = Yoshino_perovskite_conductivity(T, water_in_perovskite, P); 
    sigma_fp = Yoshino_ferropericlase_conductivity(T, water_in_ferropericlase, P, Xfe_in_ferropericlase);
    sigma_capv = Yoshino_capv_conductivity(T, water_in_capv, P);

    % Calculate effective conductivity using Hashin-Shtrikman bounds
    sigma = [sigma_pv, sigma_fp, sigma_capv];
    [sigma_upper, sigma_lower] = hashin_shtrikman(mineral_fractions, sigma);
end
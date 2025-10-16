function [sigma_upper, sigma_lower] = yoshino_mtz_lower_electrical_conductivity(T, P, ring_f, g_f, C_water, water_partition_coefficients)
% YOSHINO_MTZ_LOWER_ELECTRICAL_CONDUCTIVITY
% Compute electrical conductivity of mantle transition zone LOWER part
% using Yoshino's experimental data for ringwoodite and garnet
%
% The mantle transition zone lower part (520-660 km depth) is dominated by
% ringwoodite and garnet according to pyrolitic mantle composition.
%
% INPUTS:
%   T - Temperature in Kelvin
%   P - Pressure in Pa  
%   ring_f - Volume fraction of ringwoodite (0-1)
%   g_f - Volume fraction of garnet (0-1)
%   C_water - Total water content in wt% (typical: 0.001-0.1 wt%)
%   water_partition_coefficients - Partition coefficients [ringwoodite, garnet]
%
% REFERENCES:
% Yoshino & Katsura (2013). Electrical conductivity of mantle minerals.
% Annual Review of Earth and Planetary Sciences, 41, 605-628.
%
% Yoshino T, Manthilake G, Matsuzaki T, Katsura T. 2008a. Dry mantle transition zone 
% inferred from electrical conductivity of wadsleyite and ringwoodite. Nature 451:326–29

    % Input validation
    if T <= 0
        error('Temperature T must be positive in K, but is %.6f', T);
    end
    if P < 0
        error('Pressure P must be non-negative in Pa, but is %.6f', P);
    end
    
    % Validate mineral fractions
    mineral_fractions = [ring_f, g_f];
    if any(mineral_fractions < 0)
        error('Mineral fractions must be positive: ring_f=%.6f, g_f=%.6f', ring_f, g_f);
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
    if length(water_partition_coefficients) ~= 2
        error('water_partition_coefficients must have exactly 2 elements');
    end
    if abs(sum(water_partition_coefficients) - 1.0) > 1e-10
        error('water_partition_coefficients must sum to 1.0, but sum to %.10f', sum(water_partition_coefficients));
    end

    % Standard mantle composition: Xfe = Fe/(Fe+Mg) = 0.1
    % From Yoshino & Katsura (2013), Page 6: "for mantle peridotite (X_Fe = 0.1)"
    Xfe_in_ringwoodite = 0.1; 
    
    % Calculate individual mineral conductivities with partitioned water
    water_in_ringwoodite = C_water * water_partition_coefficients(1);
    water_in_garnet = C_water * water_partition_coefficients(2);
    
    % Use Yoshino's ringwoodite conductivity model
    sigma_ring = Yoshino_ringwoodite_conductivity(T, water_in_ringwoodite, P, Xfe_in_ringwoodite); 
    sigma_g    = Yoshino_garnet_conductivity(T, water_in_garnet, P);

    % Calculate effective conductivity using Hashin-Shtrikman bounds
    sigma = [sigma_ring, sigma_g];
    [sigma_upper, sigma_lower] = hashin_shtrikman(mineral_fractions, sigma);

end
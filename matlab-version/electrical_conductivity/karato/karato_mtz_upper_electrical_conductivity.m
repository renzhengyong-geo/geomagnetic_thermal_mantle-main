function [sigma_upper, sigma_lower] = karato_mtz_upper_electrical_conductivity(T, P, wad_f, g_f, C_water, water_partition_coefficients)
% YOSHINO_MTZ_UPPER_ELECTRICAL_CONDUCTIVITY
% Compute electrical conductivity of mantle transition zone upper part
% using Yoshino's experimental data for wadsleyite and garnet
%
% INPUTS:
%   T - Temperature in Kelvin
%   P - Pressure in Pa  
%   wad_f - Volume fraction of wadsleyite (0-1)
%   g_f - Volume fraction of garnet (0-1)
%   C_water - Total water content in wt% (typical: 0.001-0.1 wt%)
%   water_partition_coefficients - Partition coefficients [wadsleyite, garnet]
%
% REFERENCES:
% Yoshino & Katsura (2013). Electrical conductivity of mantle minerals.
% Annual Review of Earth and Planetary Sciences, 41, 605-628.

    % Input validation
    if T <= 0
        error('Temperature T must be positive, but is %.6f', T);
    end
    if P < 0
        error('Pressure P must be non-negative, but is %.6f', P);
    end
    
    % Validate mineral fractions
    mineral_fractions = [wad_f, g_f];
    if any(mineral_fractions < 0)
        error('Mineral fractions must be positive: wad_f=%.6f, g_f=%.6f', wad_f, g_f);
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

    % Calculate individual mineral conductivities with partitioned water
    water_in_wadsleyite = C_water * water_partition_coefficients(1);
    water_in_garnet = C_water * water_partition_coefficients(2);
    
    sigma_wad = Karato_wadsleyite_conductivity(T, water_in_wadsleyite, P);
    sigma_g = Karato_garnet_conductivity(T, water_in_garnet, P);

    % Calculate effective conductivity using Hashin-Shtrikman bounds
    sigma = [sigma_wad, sigma_g];
    [sigma_upper, sigma_lower] = hashin_shtrikman(mineral_fractions, sigma);

end
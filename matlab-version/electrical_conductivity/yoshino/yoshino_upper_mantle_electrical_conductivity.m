function [sigma_upper, sigma_lower]= yoshino_upper_mantle_electrical_conductivity(T, P, o_f, opx_f, cpx_f, g_f, C_water,  water_partition_coefficients)
% T is the temperature in unit of K
% P is the pressure in GPa
% o_f is the volume fraction for olivine
% opx_f is the volume fraction for opx
% cpx_f is the volume fraction for cpx
% g_f is the volume fraction for garnet
% C_water is the total water content in this mantle assemblages.
% water_partition_coefficients is the water partition coefficients
% between the olivine-opx-cpx-garnet

% Ensure T and P are positive
if T <= 0
    error('Temperature T must be positive, but is %.6f', T);
end
if P < 0
    error('Pressure P must be non-negative, but is %.6f', P);
end

% Ensure mineral fractions are positive and sum to 1.0
mineral_fractions = [o_f, opx_f, cpx_f, g_f];
if any(mineral_fractions < 0)
    error('All mineral fractions must be positive, but found negative values: o_f=%.6f, opx_f=%.6f, cpx_f=%.6f, g_f=%.6f', o_f, opx_f, cpx_f, g_f);
end
if abs(sum(mineral_fractions) - 1.0) > 1e-10
    error('Mineral fractions must sum to 1.0, but sum to %.10f: o_f=%.6f, opx_f=%.6f, cpx_f=%.6f, g_f=%.6f', sum(mineral_fractions), o_f, opx_f, cpx_f, g_f);
end

% Ensure C_water is not negative
if C_water < 0
    error('C_water must be non-negative, but is %.6f', C_water);
end

% Ensure water_partition_coefficients has exactly size 4 and sums to 1.0
if length(water_partition_coefficients) ~= 4
    error('water_partition_coefficients must have exactly 4 elements, but has %d elements', length(water_partition_coefficients));
end
if abs(sum(water_partition_coefficients) - 1.0) > 1e-10
    error('water_partition_coefficients must sum to 1.0, but sums to %.10f', sum(water_partition_coefficients));
end

% Xfe = Fe/(Fe + Mg)=0.1 in the crystal structure
% This is the standard composition used for modeling upper mantle electrical conductivity
% Yoshino, T., & Katsura, T. (2013). Electrical conductivity of mantle minerals: 
% Role of water in conductivity anomalies. Annual Review of Earth and Planetary Sciences, 
% 41, 605–628. https://doi.org/10.1146/annurev-earth-050212-124022
Xfe_in_olivine = 0.1; 
sigma_o   = Yoshino_olivine_conductivity(T, C_water*water_partition_coefficients(1), P, Xfe_in_olivine); 
sigma_opx = Yoshino_orthopyroxene_conductivity(T, C_water*water_partition_coefficients(2), P);
sigma_cpx = Yoshino_clinopyroxene_conductivity(T, C_water*water_partition_coefficients(3), P);
sigma_g   = Yoshino_garnet_conductivity(T, C_water*water_partition_coefficients(4), P);

sigma= [sigma_o, sigma_opx, sigma_cpx, sigma_g];
[sigma_upper, sigma_lower] = hashin_shtrikman(mineral_fractions, sigma);

end
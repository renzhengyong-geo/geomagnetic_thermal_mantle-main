function [T, q] = plate_cooling(depth, age, parameters)
% PLATE_COOLING - Calculate temperature and heat flux for plate cooling model
%
% MATHEMATICAL FORMULA:
%   For z ≤ L: T(z,t) = T_surface + (T_mantle - T_surface) * [z/L + (2/π) * Σ (1/n) * exp(-κn²π²t/L²) * sin(nπz/L)]
%   For z > L: T(z,t) = T_mantle (constant mantle temperature)
%
%   q(t) = k * (T_mantle - T_surface) / L * [1 + 2 * Σ exp(-κn²π²t/L²)]
%
% where:
%   L = plate thickness (m)
%   n = summation index (1 to ∞)
%   Σ = summation over n terms
%   κ = thermal diffusivity (m²/s)
%   k = thermal conductivity (W/m/K)
%
% Reference: Parsons & Sclater (1977), JGR

    % Extract parameters
    T_surface = parameters.T_surface;
    T_mantle = parameters.T_mantle;
    kappa = parameters.kappa;
    k_thermal = parameters.k_thermal;
    plate_thickness = parameters.plate_thickness;
    n_terms = parameters.n_terms;
    
    % Initialize temperature array
    T = zeros(size(depth));
    
    % Use minimum age to avoid singularity
    min_age = 1000; % 1000 seconds as minimum age
    effective_age = max(age, min_age);
       
    % Plate cooling solution with series expansion
    for i = 1:length(depth)
        z = depth(i);
        
        % If depth exceeds plate thickness, set to mantle temperature
        if z > plate_thickness
            T(i) = T_mantle;
            continue;
        end
        
        series_sum = 0;
        
        % Calculate series sum: Σ (1/n) * exp(-κn²π²t/L²) * sin(nπz/L)
        for n = 1:n_terms
            exponential_term = exp(-kappa * n^2 * pi^2 * effective_age / plate_thickness^2);
            sine_term = sin(n * pi * z / plate_thickness);
            series_sum = series_sum + (1/n) * exponential_term * sine_term;
        end
        
        % Temperature formula with series expansion (only valid for z ≤ L)
        T(i) = T_surface + (T_mantle - T_surface) * ...
               (z/plate_thickness + (2/pi) * series_sum);
    end
    
    % Calculate surface heat flux: q(t) = k(T_m - T_s)/L * [1 + 2Σ exp(-κn²π²t/L²)]
    series_sum_flux = 0;
    for n = 1:n_terms
        exponential_term = exp(-kappa * n^2 * pi^2 * effective_age / plate_thickness^2);
        series_sum_flux = series_sum_flux + exponential_term;
    end
    
    % Heat flux formula
    q = k_thermal * (T_mantle - T_surface) / plate_thickness * ...
        (1 + 2 * series_sum_flux);
end
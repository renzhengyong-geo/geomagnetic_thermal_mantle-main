function [T, q] = halfspace_cooling(depth, age, parameters)
% HALFSPACE_COOLING - Calculate temperature and heat flux for half-space cooling model
%
% MATHEMATICAL FORMULA:
%   T(z,t) = T_surface + (T_mantle - T_surface) * erf( z / (2√(κt)) )
%   q(t) = k * (T_mantle - T_surface) / √(πκt)
%
% where:
%   erf() = error function
%   z = depth (m)
%   t = lithospheric age (s)
%   κ = thermal diffusivity (m²/s)
%   k = thermal conductivity (W/m/K)
%
% Reference: Turcotte & Schubert (2014), Geodynamics

    % Extract parameters
    T_surface = parameters.T_surface;
    T_mantle = parameters.T_mantle;
    kappa = parameters.kappa;
    k_thermal = parameters.k_thermal;
    
    % Initialize temperature array
    T = zeros(size(depth));
    
    % Half-space cooling solution using error function
    for i = 1:length(depth)
        if age > 0
            % Calculate similarity variable η = z / (2√(κt))
            eta = depth(i) / (2 * sqrt(kappa * age));
            
            % Temperature formula: T(z,t) = T_s + (T_m - T_s) * erf(η)
            T(i) = T_surface + (T_mantle - T_surface) * erf(eta);
        else
            T(i) = T_surface; % Zero age case
        end
    end
    
    % Calculate surface heat flux: q(t) = k(T_m - T_s) / √(πκt)
    if age > 0
        q = k_thermal * (T_mantle - T_surface) / sqrt(pi * kappa * age);
    else
        q = inf; % Infinite heat flux at zero age
    end
end
function [sigma_HS_upper, sigma_HS_lower] = hashin_shtrikman(f, sigma, tol)
    % HASHIN_SHTRIKMAN_BOUNDS Calculate Hashin-Shtrikman bounds for effective conductivity
    %
    % References:
    % Khan, A. (2015). On Earth's mantle constitution and structure from
    %   joint analysis of geophysical and laboratory-based data: An example.
    %   Surveys in Geophysics, 36, 877-933. https://doi.org/10.1007/s10712-015-9353-z
    %
    % Original Hashin-Shtrikman theory:
    % Hashin, Z., & Shtrikman, S. (1962). A variational approach to the theory 
    %   of the effective magnetic permeability of multiphase materials. 
    %   Journal of Applied Physics, 33(10), 3125-3131.
    %
    % Mathematical formulation (Duba 1962, Khan 2015):
    % σ_HS± = [∑(x_i/(σ_i + 2σ_±))]⁻¹ - 2σ_±
    % where σ_± corresponds to min/max conductivity of the N mineral phases present
    %
    % Physical interpretation:
    % - Lower bound (HS-): non-interconnected conductive inclusions in resistive matrix  
    % - Upper bound (HS+): non-interconnected resistive inclusions in conductive matrix
    %
    % Inputs:
    %   f     - Volume fractions (vector)
    %   sigma - Component conductivities (vector)
    %   tol   - Tolerance for volume fraction sum (optional, default = 1e-6)
    %
    % Outputs:
    %   sigma_HS_upper - Upper Hashin-Shtrikman bound (HS+)
    %   sigma_HS_lower - Lower Hashin-Shtrikman bound (HS-)
    
    % Set default tolerance
    if nargin < 3
        tol = 1e-6;
    end
    
    % Validate inputs
    assert(all(f >= 0), 'Error: Volume fractions must be non-negative.');
    assert(abs(sum(f) - 1) < tol, 'Error: Volume fractions must sum to 1 within tolerance.');
    assert(all(sigma > 0), 'Error: Conductivities must be positive.');
    assert(length(f) == length(sigma), 'Error: f and sigma must have same length.');
    
    n_phases = length(f);
    
    % Handle single-phase material
    if n_phases == 1
        sigma_HS_upper = sigma(1);
        sigma_HS_lower = sigma(1);
        return;
    end
    
    % Identify max and min conductivities for σ_±
    sigma_plus = max(sigma);   % σ_+ for upper bound
    sigma_minus = min(sigma);  % σ_- for lower bound
    
    % Compute upper bound (HS+)
    % σ_HS+ = [∑(x_i/(σ_i + 2σ_+))]⁻¹ - 2σ_+
    sum_upper = 0;
    for i = 1:n_phases
        sum_upper = sum_upper + f(i) / (sigma(i) + 2 * sigma_plus);
    end
    
    sigma_HS_upper = 1/sum_upper - 2*sigma_plus;
    
    % Compute lower bound (HS-)
    % σ_HS- = [∑(x_i/(σ_i + 2σ_-))]⁻¹ - 2σ_-
    sum_lower = 0;
    for i = 1:n_phases
        sum_lower = sum_lower + f(i) / (sigma(i) + 2 * sigma_minus);
    end
    
    sigma_HS_lower = 1/sum_lower - 2*sigma_minus;
    
    % Ensure bounds are physically reasonable
    sigma_HS_upper = max(sigma_HS_upper, sigma_HS_lower);
end


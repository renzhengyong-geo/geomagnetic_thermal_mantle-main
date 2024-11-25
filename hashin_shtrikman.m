function [sigma_HS_upper, sigma_HS_lower] = hashin_shtrikman(f, sigma)
    % Validate inputs
    assert(all(f >= 0), 'Error: Volume fractions must be non-negative.');
    % assert(abs(sum(f) - 1) < 1e-6, 'Error: Volume fractions must sum to 1.');
    assert(all(sigma > 0), 'Error: Conductivities must be positive.');

    % Identify max and min conductivities
    [sigma_n, idx_max] = max(sigma);
    [sigma_1, idx_min] = min(sigma);

    % Compute A_plus for upper bound
    A_plus = 0;
    for i = 1:length(sigma)
        if i ~= idx_max
            diff = sigma(i) - sigma_n;
            assert(diff < 0, 'Error: Conductivity difference (sigma_i - sigma_n) must be negtive.');
            A_plus = A_plus + f(i) / ((1 / diff) + (1 / (3 * sigma_n)));
        end
    end
    denominator_upper = 1 - A_plus / (3 * sigma_n);
    assert(denominator_upper > 0, 'Error: Denominator for HS upper bound is invalid.');
    sigma_HS_upper = sigma_n + A_plus / denominator_upper;

    % Compute A_minus for lower bound
    A_minus = 0;
    for i = 1:length(sigma)
        if i ~= idx_min
            diff = sigma(i) - sigma_1;
            assert(diff > 0, 'Error: Conductivity difference (sigma_i - sigma_1) must be positive.');
            A_minus = A_minus + f(i) / ((1 / diff) + (1 / (3 * sigma_1)));
        end
    end
    denominator_lower = 1 - A_minus / (3 * sigma_1);
    assert(denominator_lower > 0, 'Error: Denominator for HS lower bound is invalid.');
    sigma_HS_lower = sigma_1 + A_minus / denominator_lower;
end

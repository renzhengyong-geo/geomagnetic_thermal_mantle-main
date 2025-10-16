function [P, T, aggregate_sigma_plus, aggregate_sigma_minus] = ...
    compute_lower_mantle_conductivity_depth_profile(lowerMantleTable, water_content, group_id)
    % Compute conductivity depth profile for the lower mantle.
    %
    % INPUTS:
    %   lowerMantleTable         - Table with columns: 'P', 'T', 'Pv', 'Ca_pv', 'Fp', 'Ppv'
    %   water_content            - Water content for Pv (Perovskite) (wt%).
    %   group_id                 - Group ('Karato' or 'Yashino') for conductivity calculation.
    %
    % OUTPUTS:
    %   P                    - Pressure profile (bar).
    %   T                    - Temperature profile (K).
    %   aggregate_sigma_plus - Upper Hashin-Shtrikman bound (S/m).
    %   aggregate_sigma_minus- Lower Hashin-Shtrikman bound (S/m).

    % -------------------------------------------------------------------------
    % Input Validation
    % -------------------------------------------------------------------------
    m = size(lowerMantleTable, 1); % Number of rows

    % Check table structure
    requiredColumns = {'P', 'T', 'Pv', 'Ca_pv', 'Fp', 'Ppv'};
    missingColumns = setdiff(requiredColumns, lowerMantleTable.Properties.VariableNames);
    assert(isempty(missingColumns), 'Error: Missing required columns: %s.', strjoin(missingColumns, ', '));
    assert(m == length(water_content), 'Error: Rows in lowerMantleTable must match length of water_content.');
    assert(ismember(lower(group_id), {'karato', 'yashino'}), 'Error: group_id must be "karato" or "yashino".');

    % -------------------------------------------------------------------------
    % Extract Data
    % -------------------------------------------------------------------------
    P = lowerMantleTable.P;       % Pressure (bar)
    T = lowerMantleTable.T;       % Temperature (K)
    pv = lowerMantleTable.Pv;     % Perovskite volume fraction
    capv = lowerMantleTable.Ca_pv; % Calcium perovskite volume fraction
    fp = lowerMantleTable.Fp;     % Ferropericlase volume fraction
    ppv = lowerMantleTable.Ppv;   % Post-perovskite volume fraction

    % Volume fractions matrix
    f = [pv, capv, fp, ppv];

    % Validate volume fractions
    assert(all(sum(f, 2) > 0), 'Error: Total volume fractions must be positive for each row.');

    % Normalize volume fractions
    f = f ./ sum(f, 2);

    % -------------------------------------------------------------------------
    % Calculate Water Partitioning
    % -------------------------------------------------------------------------
    D = ones(m, 3); % Partitioning coefficients for Ca-Pv, Fp, and PPV

    % Calculate water content for each mineral
    Cw = zeros(m, 4); % Water content for Pv, Ca-Pv, Fp, and PPV
    Cw(:, 1) = water_content; % Water content for Pv
    Cw(:, 2) = Cw(:, 1) ./ D(:, 1); % Water content for Ca-Pv
    Cw(:, 3) = Cw(:, 1) ./ D(:, 2); % Water content for Fp
    Cw(:, 4) = Cw(:, 1) ./ D(:, 3); % Water content for PPV

    % Validate water content
    assert(all(Cw(:) >= 0), 'Error: Water content (Cw) must be non-negative.');

    % -------------------------------------------------------------------------
    % Compute Conductivities
    % -------------------------------------------------------------------------
    sigma = zeros(m, 4); % Conductivities for Pv, Ca-Pv, Fp, and PPV
    for i = 1:m
        t = T(i); % Temperature (K)
        p = P(i) * 1e5; % Pressure (Pa)
        switch lower(group_id)
            case 'karato'
                sigma(i, 1) = Karato_perovskite_conductivity(t, Cw(i, 1), p);
                sigma(i, 2) = Karato_calcium_perovskite_conductivity(t, Cw(i, 2), p);
                sigma(i, 3) = Karato_ferropericlase_conductivity(t, Cw(i, 3), p);
                sigma(i, 4) = Karato_post_perovskite_conductivity(t, Cw(i, 4), p);
            case 'yashino'
                sigma(i, 1) = Yashino_perovskite_conductivity(t, Cw(i, 1), p);
                sigma(i, 2) = Yashino_calcium_perovskite_conductivity(t, Cw(i, 2), p);
                sigma(i, 3) = Yashino_ferropericlase_conductivity(t, Cw(i, 3), p);
                sigma(i, 4) = Yashino_post_perovskite_conductivity(t, Cw(i, 4), p);
            otherwise
                error('Unexpected group_id value: %s', group_id);
        end
    end

    % Validate conductivities
    assert(all(sigma(:) > 0), 'Error: Conductivities must be positive.');

    % -------------------------------------------------------------------------
    % Compute Aggregate Conductivities
    % -------------------------------------------------------------------------
    aggregate_sigma_plus = zeros(m, 1);
    aggregate_sigma_minus = zeros(m, 1);

    for i = 1:m
        fractions = f(i, :); % Volume fractions
        conductivities = sigma(i, :); % Conductivities (S/m)
        % Ensure valid inputs for Hashin-Shtrikman
        [sigma_HS_upper, sigma_HS_lower] = hashin_shtrikman(fractions, conductivities);
        aggregate_sigma_plus(i) = sigma_HS_upper;
        aggregate_sigma_minus(i) = sigma_HS_lower;
    end
end

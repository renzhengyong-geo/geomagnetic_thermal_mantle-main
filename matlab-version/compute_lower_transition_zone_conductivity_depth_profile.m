function [P, T, aggregate_sigma_plus, aggregate_sigma_minus,total_Cw] = ...
    compute_lower_transition_zone_conductivity_depth_profile(lowerTransitionZoneTable, water_content, group_id)
    % Compute conductivity depth profile for the lower transition zone.
    %
    % INPUTS:
    %   lowerTransitionZoneTable - Table with columns: P, T, Ring, Gt
    %   water_content            - Water content for Ring(Ringwoodite) (wt%).
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
    m = size(lowerTransitionZoneTable, 1); % Number of rows

    % Check table structure
    requiredColumns = {'P', 'T', 'Ring', 'Gt'};
    missingColumns = setdiff(requiredColumns, lowerTransitionZoneTable.Properties.VariableNames);
    assert(isempty(missingColumns), 'Error: Missing required columns: %s.', strjoin(missingColumns, ', '));
    assert(m == length(water_content), 'Error: Rows in lowerTransitionZoneTable must match length of water_content.');
    assert(ismember(lower(group_id), {'karato', 'yashino'}), 'Error: group_id must be "karato" or "yashino".');

    % -------------------------------------------------------------------------
    % Extract Data
    % -------------------------------------------------------------------------
    P = lowerTransitionZoneTable.P;       % Pressure (bar)
    T = lowerTransitionZoneTable.T;       % Temperature (K)
    ring = lowerTransitionZoneTable.Ring; % Ringwoodite volume fraction
    gt = lowerTransitionZoneTable.Gt;     % Garnet volume fraction

    % Volume fractions matrix
    f = [ring, gt];

    % Validate volume fractions
    assert(all(sum(f, 2) > 0), 'Error: Volume fractions must sum to a positive value for each row.');

    % -------------------------------------------------------------------------
    % Calculate Water Partitioning
    %   Based on: Inoue, T., Katsuda, M., Yurimoto, H., & Irifune, T. (2006).
    % -------------------------------------------------------------------------
    D = ones(m, 1) * 5; % Partitioning for gt (hardcoded as 5 for this example)

    % Calculate water content for each mineral
    Cw = zeros(m, 2); % Water content for Ring, Gt
    Cw(:, 1) = water_content; % Water content for Ring
    Cw(:, 2) = Cw(:, 1) ./ D; % Distribute water content to Gt based on partitioning

    % Validate water content
    assert(all(Cw(:) >= 0), 'Error: Water content (Cw) must be non-negative.');

    % -------------------------------------------------------------------------
    % Compute Conductivities
    % -------------------------------------------------------------------------
    sigma = zeros(m, 2); % Conductivities for Ring, Gt
    for i = 1:m
        t = T(i); % Temperature (K)
        p = P(i) * 1e5; % Pressure (Pa)
        switch lower(group_id)
            case 'karato'
                sigma(i, 1) = Karato_ringwoodite_conductivity(t, Cw(i, 1), p);
                sigma(i, 2) = Karato_garnet_conductivity(t, Cw(i, 2), p);
            case 'yashino'
                sigma(i, 1) = Yashino_ringwoodite_conductivity(t, Cw(i, 1), p);
                sigma(i, 2) = Yashino_garnet_conductivity(t, Cw(i, 2), p);
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
        fractions = f(i, :) / 100; % Convert to fractions
        conductivities = sigma(i, :); % Conductivities (S/m)
        [sigma_HS_upper, sigma_HS_lower] = hashin_shtrikman(fractions, conductivities);
        aggregate_sigma_plus(i) = sigma_HS_upper;
        aggregate_sigma_minus(i) = sigma_HS_lower;
    end

    % -------------------------------------------------------------------------
    % Compute Aggregate Water content
    % -------------------------------------------------------------------------
    total_Cw = zeros(m, 1);
    for i = 1:m
        for j=1:2
          total_Cw(i) = total_Cw(i)+ (f(i, j) / 100)*Cw(i, j); % Convert to fractions
        end
    end
end

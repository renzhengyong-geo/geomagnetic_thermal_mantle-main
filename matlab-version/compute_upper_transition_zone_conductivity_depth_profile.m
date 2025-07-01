function [P, T, aggregate_sigma_plus, aggregate_sigma_minus,total_Cw] = ...
    compute_upper_transition_zone_conductivity_depth_profile(upperTransitionZoneTable, water_content, group_id)
    % Compute conductivity depth profile for upper_transition_zone
    %
    % INPUTS:
    %   upperTransitionZoneTable - Table with columns: P, T, Wad, Gt, Ring
    %   water_content            - Water content for Wad(wadsleyite) (wt%).
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
    m = size(upperTransitionZoneTable, 1); % Number of rows
    % Check table structure
    requiredColumns = {'P', 'T', 'Wad', 'Gt', 'Ring'};
    missingColumns = setdiff(requiredColumns, upperTransitionZoneTable.Properties.VariableNames);
    assert(isempty(missingColumns), 'Error: Missing required columns: %s.', strjoin(missingColumns, ', '));
    assert(m == length(water_content), 'Error: Rows in upperTransitionZoneTable must match length of water_content.');
    assert(ismember(lower(group_id), {'karato', 'yoshino'}), 'Error: group_id must be "karato" or "yashino".');

    % -------------------------------------------------------------------------
    % Extract Data
    % -------------------------------------------------------------------------
    P = upperTransitionZoneTable.P;       % Pressure (bar)
    T = upperTransitionZoneTable.T;       % Temperature (K)
    wad = upperTransitionZoneTable.Wad;   % Wadsleyite volume fraction
    gt = upperTransitionZoneTable.Gt;     % Garnet volume fraction
    ring = upperTransitionZoneTable.Ring; % Ringwoodite volume fraction

    % Volume fractions matrix
    f = [wad, gt, ring];

    % Validate volume fractions
    assert(all(sum(f, 2) > 0), 'Error: Volume fractions must sum to a positive value for each row.');

    % -------------------------------------------------------------------------
    % Calculate Water Partitioning
    % -------------------------------------------------------------------------
    D = zeros(m, 2); % Partitioning for gt and ring
    for i = 1:m
        output = water_partition_wad_ring_gt(T(i) - 273.15, P(i) / 1000); % Convert K to °C, bar to kbar
        D(i, :) = [output(4), output(3)]; % Partitioning coefficients: Wad vs Gt, Wad vs Ring
    end

    % Calculate water content for each mineral
    Cw = zeros(m, 3); % Water content for Wad, Gt, Ring
    Cw(:, 1) = water_content; % Water content for Wad
    for i = 2:3
        Cw(:, i) = Cw(:, 1) ./ D(:, i - 1); % Distribute water content to Gt and Ring based on partitioning
    end

    % -------------------------------------------------------------------------
    % Compute Conductivities
    % -------------------------------------------------------------------------
    sigma = zeros(m, 3); % Conductivities for Wad, Gt, and Ring 
    for i = 1:m
        t = T(i); % Temperature (K)
        p = P(i) * 1e5; % Pressure (Pa)
        switch lower(group_id)
            case 'karato'
                sigma(i, 1) = Karato_wadsleyite_conductivity(t, Cw(i, 1), p);
                sigma(i, 2) = Karato_garnet_conductivity(t, Cw(i, 2), p);
                sigma(i, 3) = Karato_ringwoodite_conductivity(t, Cw(i, 3), p);
            case 'yoshino'
                sigma(i, 1) = Yoshino_wadsleyite_conductivity(t, Cw(i, 1), p, 1.0);
                sigma(i, 2) = Yoshino_garnet_conductivity(t, Cw(i, 2), p);
                sigma(i, 3) = Yoshino_ringwoodite_conductivity(t, Cw(i, 3), p, 1.0);
            otherwise
                error('Unexpected group_id value: %s', group_id);
        end
    end

    % -------------------------------------------------------------------------
    % Compute Aggregate Conductivities
    % -------------------------------------------------------------------------
    aggregate_sigma_plus = zeros(m, 1);
    aggregate_sigma_minus = zeros(m, 1);
    
    Cw
    f
    sigma
    
    for i = 1:m
        fractions = f(i, :) / 100; % Convert to fractions
        conductivities = sigma(i, :); % Conductivities (S/m)
        % Ensure valid inputs for Hashin-Shtrikman
        assert(all(conductivities > 0), 'Error: Conductivities must be positive.');
        % assert(abs(sum(fractions) - 1) < 1e-6, 'Error: Volume fractions must sum to 1.');
        [sigma_HS_upper, sigma_HS_lower] = hashin_shtrikman(fractions, conductivities);
        aggregate_sigma_plus(i) = sigma_HS_upper;
        aggregate_sigma_minus(i) = sigma_HS_lower;
    end

    % -------------------------------------------------------------------------
    % Compute Aggregate Water content
    % -------------------------------------------------------------------------
    total_Cw = zeros(m, 1);
    for i = 1:m
        for j=1:3
          total_Cw(i) = total_Cw(i)+ (f(i, j) / 100)*Cw(i, j); % Convert to fractions
        end
    end

end

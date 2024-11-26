function [P, T, aggregate_sigma_plus, aggregate_sigma_minus, total_Cw] = compute_upper_mantle_conductivity_depth_profile(upperMantleTable, ...
    water_content, group_id)
    % Compute conductivity depth profile for the upper mantle.
    %
    % INPUTS:
    %   upperMantleTable - Table with columns: P, T, O, Cpx, Opx, Gt.
    %   water_content    - Water content for O(olivine) (wt%).
    %   group_id         - Group ('Karato' or 'Yashino') for conductivity calculation.
    %
    % OUTPUTS:
    %   P                    - Pressure profile (bar).
    %   T                    - Temperature profile (K).
    %   aggregate_sigma_plus - Upper Hashin-Shtrikman bound (S/m).
    %   aggregate_sigma_minus- Lower Hashin-Shtrikman bound (S/m).

    % -------------------------------------------------------------------------
    % Input Validation
    % -------------------------------------------------------------------------
    m = size(upperMantleTable, 1); % Number of rows
    assert(size(upperMantleTable, 2) == 6, 'Error: Columns in upperMantleTable must be 6 (P, T, O, Cpx, Opx, Gt).');
    assert(m == length(water_content), 'Error: Rows in upperMantleTable must match length of water_content.');
    assert(ismember(lower(group_id), {'karato', 'yashino'}), 'Error: group_id must be "karato" or "yashino".');

    % -------------------------------------------------------------------------
    % Extract Data
    % -------------------------------------------------------------------------
    P       = upperMantleTable.P;      % Pressure (bar)
    T       = upperMantleTable.T;      % Temperature (K)
    olivine = upperMantleTable.O;      % Olivine volume fraction
    cpx     = upperMantleTable.Cpx;    % Clinopyroxene volume fraction
    opx     = upperMantleTable.Opx;    % Orthopyroxene volume fraction
    gt      = upperMantleTable.Gt;     % Garnet volume fraction

    % Volume fractions matrix
    f = [olivine, cpx, opx, gt];

    % -------------------------------------------------------------------------
    % Calculate Water Partitioning
    % -------------------------------------------------------------------------
    D = zeros(m, 3); % Partitioning for Cpx, Opx, Gt
    for i = 1:m
        output = water_partition(T(i) - 273.15, P(i) / 1000); % Convert K to °C, bar to kbar
        D(i, :) = output; % Partitioning coefficients
    end

    % Calculate water content for each mineral
    Cw = zeros(m, 4); % Water content for O, Cpx, Opx, Gt
    Cw(:, 1) = water_content; % Olivine
    for i = 2:4
        Cw(:, i) = Cw(:, 1) ./ D(:, i - 1); % Cpx, Opx, Gt
    end

    % -------------------------------------------------------------------------
    % Compute Conductivities
    % -------------------------------------------------------------------------
    sigma = zeros(m, 4); % Conductivities for O, Cpx, Opx, Gt
    for i = 1:m
        t = T(i); % Temperature (K)
        p = P(i) * 1e5; % Pressure (Pa)
        switch lower(group_id)
            case 'karato'
                sigma(i, 1) = Karato_olivine_conductivity(t, Cw(i, 1), p);
                sigma(i, 2) = Karato_clinopyroxene_conductivity(t, Cw(i, 2), p);
                sigma(i, 3) = Karato_orthopyroxene_conductivity(t, Cw(i, 3), p);
                sigma(i, 4) = Karato_garnet_conductivity(t, Cw(i, 4), p);
            case 'yashino'
                % sigma(i, 1) = Yashino_olivine_conductivity(t, Cw(i, 1), p);
                % sigma(i, 2) = Yashino_clinopyroxene_conductivity(t, Cw(i, 2), p);
                % sigma(i, 3) = Yashino_orthopyroxene_conductivity(t, Cw(i, 3), p);
                % sigma(i, 4) = Yashino_garnet_conductivity(t, Cw(i, 4), p);
            otherwise
                error('Unexpected group_id value: %s', group_id);
        end
    end

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
        for j=1:4
          total_Cw(i) = total_Cw(i)+ (f(i, j) / 100)*Cw(i, j); % Convert to fractions
        end
    end

end

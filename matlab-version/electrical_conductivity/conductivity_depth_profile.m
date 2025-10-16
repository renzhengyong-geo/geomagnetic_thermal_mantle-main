function [Pressure, Temperature, Sigma_plus, Sigma_minus] = conductivity_depth_profile(uppermantle_layers, transiztionzone_layers, lowermantle_layers, columnNames, dataMatrix, water_content, group_id)
    % conductivity_depth_profile Computes the depth and conductivity profile.
    %
    % INPUTS:
    %   uppermantle_layers      - Number of layers in the upper mantle.
    %   transiztionzone_layers  - Number of layers in the transition zone.
    %   lowermantle_layers      - Number of layers in the lower mantle.
    %   columnNames             - Cell array of column names from the header.
    %   dataMatrix              - Numeric matrix containing the data.
    %   water_content           - Water content for relevant minerals (wt%).
    %   group_id                - Group name ('Karato' or 'Yashino') to select the routine.
    %
    % OUTPUTS:
    %   Pressure    - Pressure profile for the upper mantle.
    %   Sigma_plus  - Upper Hashin-Shtrikman conductivity profile.
    %   Sigma_minus - Lower Hashin-Shtrikman conductivity profile.
    %
    % NOTES:
    %   The code assumes that the input data matrix contains columns labeled with
    %   the required mineral and phase data.

    % -------------------------------------------------------------------------
    % Input Validation
    % -------------------------------------------------------------------------
    % Check that the number of layers matches the number of rows in dataMatrix
    n_layers = uppermantle_layers + transiztionzone_layers + lowermantle_layers;
    assert(n_layers == size(dataMatrix, 1), ...
        'Error: n_layers does not match the number of rows in dataMatrix.');

    % Check that the water_content vector matches the number of layers
    assert(n_layers == length(water_content), ...
        'Error: n_layers does not match the length of water_content.');

    % Validate group_id
    assert(ismember(lower(group_id), {'karato', 'yashino'}), ...
        'Error: group_id must be "karato" or "yashino" (case-insensitive).');

    % Ensure required columns are present in columnNames
    pressureIdx = find(strcmpi(columnNames, 'P(bar)'), 1);
    temperatureIdx = find(strcmpi(columnNames, 'T(K)'), 1);
    if isempty(pressureIdx) || isempty(temperatureIdx)
        error('Missing required columns in columnNames. Ensure "P(bar)" and "T(K)" are present.');
    end

    % -------------------------------------------------------------------------
    % Extract Data for Upper Mantle
    % -------------------------------------------------------------------------
    % Extract temperature and pressure
    pressure = dataMatrix(:, pressureIdx);       % Pressure in bar
    temperature = dataMatrix(:, temperatureIdx); % Temperature in Kelvin

    % Extract volume fractions for upper mantle minerals
    upper_mantle_minerals = {'O', 'Cpx', 'Opx', 'Gt'};
    upper_mantle_volumeFractions = zeros(uppermantle_layers, numel(upper_mantle_minerals));
    for i = 1:numel(upper_mantle_minerals)
        idx = find(strcmpi(columnNames, upper_mantle_minerals{i}), 1);
        if ~isempty(idx)
            upper_mantle_volumeFractions(:, i) = dataMatrix(1:uppermantle_layers, idx);
        else
            warning('Mineral "%s" not found in columnNames.', upper_mantle_minerals{i});
        end
    end

    % Extract water content and calculate partitioning
    upper_mantle_water_partition_table = zeros(uppermantle_layers, 3); % For Cpx, Opx, Gt
    upper_mantle_water_content = zeros(uppermantle_layers, 4);         % For O, Cpx, Opx, Gt
    olivine_water_content = water_content(1:uppermantle_layers);

    for i = 1:uppermantle_layers
        output = water_partition(temperature(i) - 273.15, pressure(i) / 1000); % Convert K to °C, bar to kbar
        upper_mantle_water_partition_table(i, :) = output; % Partitioning for Cpx, Opx, Gt
        upper_mantle_water_content(i, 1) = olivine_water_content(i); % O
        upper_mantle_water_content(i, 2) = olivine_water_content(i) / output(2); % Cpx
        upper_mantle_water_content(i, 3) = olivine_water_content(i) / output(1); % Opx
        upper_mantle_water_content(i, 4) = olivine_water_content(i) / output(3); % Gt
    end

    % -------------------------------------------------------------------------
    % Compute Conductivities for Upper Mantle Minerals
    % -------------------------------------------------------------------------
    upper_mantle_conductivity = zeros(uppermantle_layers, 4); % For O, Cpx, Opx, Gt
    for i = 1:uppermantle_layers
        T = temperature(i);  % Temperature in Kelvin
        P = pressure(i) * 1e5; % Convert bar to Pascal

        % Select the correct conductivity function based on group_id
        switch lower(group_id)
            case 'karato'
                upper_mantle_conductivity(i, 1) = Karato_olivine_conductivity(T, upper_mantle_water_content(i, 1), P);
                upper_mantle_conductivity(i, 2) = Karato_clinopyroxene_conductivity(T, upper_mantle_water_content(i, 2), P);
                upper_mantle_conductivity(i, 3) = Karato_orthopyroxene_conductivity(T, upper_mantle_water_content(i, 3), P);
                upper_mantle_conductivity(i, 4) = Karato_garnet_conductivity(T, upper_mantle_water_content(i, 4), P);
            case 'yashino'
                % upper_mantle_conductivity(i, 1) = Yashino_olivine_conductivity(T, upper_mantle_water_content(i, 1), P);
                % upper_mantle_conductivity(i, 2) = Yashino_clinopyroxene_conductivity(T, upper_mantle_water_content(i, 2), P);
                % upper_mantle_conductivity(i, 3) = Yashino_orthopyroxene_conductivity(T, upper_mantle_water_content(i, 3), P);
                % upper_mantle_conductivity(i, 4) = Yashino_garnet_conductivity(T, upper_mantle_water_content(i, 4), P);
            otherwise
                error('Unexpected group_id value: %s', group_id);
        end
    end
    % -------------------------------------------------------------------------
    % Compute Aggregate Conductivity Using Hashin-Shtrikman Bounds
    % -------------------------------------------------------------------------
    upper_mantle_aggregate_sigma_plus = zeros(uppermantle_layers, 1);
    upper_mantle_aggregate_sigma_minus = zeros(uppermantle_layers, 1);
    disp("upper_mantle_conductivity:"); 
    disp(upper_mantle_conductivity);
    disp("upper_mantle_volumeFractions:"); 
    disp(upper_mantle_volumeFractions);
    for i = 1:uppermantle_layers
        f = upper_mantle_volumeFractions(i, :)/100; % Volume fractions, in %
        % % Normalize f; do a little changes to volume fractions so that
        % % satisfy the requirment.
        % f = f / sum(f);
        sigma = upper_mantle_conductivity(i, :); % Conductivities, S/m
        [sigma_HS_upper, sigma_HS_lower] = hashin_shtrikman(f, sigma);
        upper_mantle_aggregate_sigma_plus(i) = sigma_HS_upper;
        upper_mantle_aggregate_sigma_minus(i) = sigma_HS_lower;
    end

    % Assign outputs
    Pressure = pressure(1:uppermantle_layers);
    Temperature=temperature(1:uppermantle_layers);
    Sigma_plus = upper_mantle_aggregate_sigma_plus;
    Sigma_minus = upper_mantle_aggregate_sigma_minus;
end

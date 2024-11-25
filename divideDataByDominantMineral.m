function [upperMantleTable, upperTransitionZoneTable, lowerTransitionZoneTable, lowerMantleTable] = divideDataByDominantMineral(inputData, columnNames)
    % Function to divide the input data into four tables based on dominant minerals.
    %
    % INPUTS:
    %   inputData    - Input data matrix
    %   columnNames  - Cell array of column names
    %
    % OUTPUTS:
    %   upperMantleTable         - Table for the upper mantle
    %   upperTransitionZoneTable - Table for the upper part of the transition zone
    %   lowerTransitionZoneTable - Table for the lower part of the transition zone
    %   lowerMantleTable         - Table for the lower mantle

    % Identify column indices
    P_idx = find(strcmpi(columnNames, 'P(bar)'));
    T_idx = find(strcmpi(columnNames, 'T(K)'));
    O_idx = find(strcmpi(columnNames, 'O'));
    Cpx_idx = strcmpi(columnNames, 'Cpx');
    Opx_idx = strcmpi(columnNames, 'Opx');
    Gt_idx = find(strcmpi(columnNames, 'Gt'));
    Wad_idx = find(strcmpi(columnNames, 'Wad'));
    Ring_idx = find(strcmpi(columnNames, 'Ring'));
    Pv_idx = find(strcmpi(columnNames, 'Pv'));
    Ca_pv_idx = find(strcmpi(columnNames, 'Ca-pv'));
    Fp_idx = find(strcmpi(columnNames, 'Fp'));
    Ppv_idx = find(strcmpi(columnNames, 'Ppv'));

    % Create default columns with zeros if missing
    numRows = size(inputData, 1);
    defaultCol = zeros(numRows, 1);
    Gt_col = defaultCol; if ~isempty(Gt_idx), Gt_col = inputData(:, Gt_idx); end
    Wad_col = defaultCol; if ~isempty(Wad_idx), Wad_col = inputData(:, Wad_idx); end
    Ring_col = defaultCol; if ~isempty(Ring_idx), Ring_col = inputData(:, Ring_idx); end
    Ca_pv_col = defaultCol; if ~isempty(Ca_pv_idx), Ca_pv_col = inputData(:, Ca_pv_idx); end
    Fp_col = defaultCol; if ~isempty(Fp_idx), Fp_col = inputData(:, Fp_idx); end
    Ppv_col = defaultCol; if ~isempty(Ppv_idx), Ppv_col = inputData(:, Ppv_idx); end

    % Partition data into subsets
    upperMantleRows = inputData(:, O_idx) > 0;
    upperTransitionZoneRows = Wad_col > 0 & ~upperMantleRows;
    lowerTransitionZoneRows = Ring_col > 0 & ~(upperMantleRows | upperTransitionZoneRows);
    lowerMantleRows = Pv_idx > 0 & ~(upperMantleRows | upperTransitionZoneRows | lowerTransitionZoneRows);

    % Upper Mantle Table
    upperMantleData = [
        inputData(upperMantleRows, P_idx), ...
        inputData(upperMantleRows, T_idx), ...
        inputData(upperMantleRows, O_idx), ...
        inputData(upperMantleRows, Cpx_idx), ...
        inputData(upperMantleRows, Opx_idx), ...
        Gt_col(upperMantleRows)
    ];
    upperMantleTable = createTable(upperMantleData, {'P', 'T', 'O', 'Cpx', 'Opx', 'Gt'});

    % Upper Transition Zone Table
    upperTransitionZoneData = [
        inputData(upperTransitionZoneRows, P_idx), ...
        inputData(upperTransitionZoneRows, T_idx), ...
        Wad_col(upperTransitionZoneRows), ...
        Gt_col(upperTransitionZoneRows), ...
        Ring_col(upperTransitionZoneRows)
    ];
    upperTransitionZoneTable = createTable(upperTransitionZoneData, {'P', 'T', 'Wad', 'Gt', 'Ring'});

    % Lower Transition Zone Table
    lowerTransitionZoneData = [
        inputData(lowerTransitionZoneRows, P_idx), ...
        inputData(lowerTransitionZoneRows, T_idx), ...
        Ring_col(lowerTransitionZoneRows), ...
        Gt_col(lowerTransitionZoneRows)
    ];
    lowerTransitionZoneTable = createTable(lowerTransitionZoneData, {'P', 'T', 'Ring', 'Gt'});

    % Lower Mantle Table
    lowerMantleData = [
        inputData(lowerMantleRows, P_idx), ...
        inputData(lowerMantleRows, T_idx), ...
        inputData(lowerMantleRows, Pv_idx), ...
        Ca_pv_col(lowerMantleRows), ...
        Fp_col(lowerMantleRows), ...
        Ppv_col(lowerMantleRows)
    ];
    lowerMantleTable = createTable(lowerMantleData, {'P', 'T', 'Pv', 'Ca_pv', 'Fp', 'Ppv'});
end

function tableOut = createTable(data, variableNames)
    % Helper function to create a table with specified variable names
    if isempty(data)
        % Create an empty table with specified column names
        tableOut = array2table(zeros(0, numel(variableNames)), 'VariableNames', variableNames);
    else
        % Check consistency between data and variable names
        if size(data, 2) ~= numel(variableNames)
            error('Mismatch between data columns and variable names.');
        end
        % Create table with actual data
        tableOut = array2table(data, 'VariableNames', variableNames);
    end
end

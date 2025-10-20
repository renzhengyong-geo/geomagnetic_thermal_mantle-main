function fO2 = compute_fO2_Righter(bufferName, temperature_C)
    % COMPUTE_FO2 calculates the absolute fO2 for a specified buffer and temperature.
    % 
    % 1. Barin, I. (1995). Thermochemical Data of Pure Substances, 3rd ed. 
    %    Wiley-VCH Verlag GmbH.
    % 
    % 2. Righter, K., et al. (2023). Oxygen fugacity buffering in high-pressure solid 
    %    media assemblies from IW-6.5 to IW+4.5 and application to the V K-edge oxybarometer.
    %    American Mineralogist, 108(3): 498–513. 
    %    https://doi.org/10.2138/am-2022-8301

    % Inputs:
    %   bufferName - Name of the buffer (case insensitive), e.g., 'Re_ReO2'
    %   temperature_C - Temperature in degrees Celsius (°C)
    %
    % Output:
    %   fO2 - The absolute fO2 value in bar

    % Buffer constants (a and b values)
    buffers = struct( ...
        'Re_ReO2', [-22552, 9.2559], ...
        'Ni_NiO', [-24522, 8.9248], ...
        'Co_CoO', [-24445, 7.3503], ...
        'W_WO3', [-28629, 8.844], ...
        'Fe_FeO', [-27847, 6.2606], ... % IW buffer
        'Mo_MoO2', [-29894, 8.6839], ...
        'Cr_Cr2O3', [-38907, 8.4065], ...
        'V_V2O3', [-41571, 8.0921], ...
        'Ta_Ta2O5', [-41941, 8.4072], ...
        'Nb_NbO', [-43165, 8.8919], ...
        'Si_SiO2', [-47283, 9.1191]); % New curve added

    % Convert buffer name to case-insensitive match
    fieldNames = fieldnames(buffers);
    bufferName = lower(bufferName); % Convert input to lowercase
    matchedField = fieldNames(strcmpi(fieldNames, bufferName)); % Match field

    % Validate buffer name
    if isempty(matchedField)
        error('Invalid buffer name. Please provide a valid buffer name.');
    end

    % Extract buffer constants
    constants = buffers.(matchedField{1});
    a = constants(1);
    b = constants(2);

    % Convert temperature to Kelvin
    temperature_K = temperature_C + 273.15;

    % Calculate log10(fO2)
    log_fO2 = a / temperature_K + b;

    % Convert to absolute fO2
    fO2 = 10^log_fO2; % fO2 in bar
end

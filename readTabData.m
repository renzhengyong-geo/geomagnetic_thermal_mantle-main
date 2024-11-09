function [columnNames, dataMatrix] = readTabData(fileName)
    % Function to read data from a tab file.
    % Inputs:
    %   fileName - The name of the input file (string)
    % Outputs:
    %   columnNames - A string array containing column headers
    %   dataMatrix - A matrix containing numerical data with NaNs replaced by zeros

    % Open the file
    fileID = fopen(fileName, 'r');
    
    % Check if file opened successfully
    if fileID == -1
        error('Could not open the file.');
    end
    
    % Skip the initial lines (version, file name, irrelevant data)
    for i = 1:8
        fgetl(fileID);
    end

    % Read the column headers (line 9)
    headerLine = fgetl(fileID);
    columnNames = strsplit(strtrim(headerLine));
    columnNames = string(columnNames);
    
    % Read the remaining data starting from line 10 onwards
    data = textscan(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f', 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
    
    % Close the file
    fclose(fileID);
    
    % Convert the cell array to a matrix
    dataMatrix = cell2mat(data);
    
    % Replace NaN values with zeros
    dataMatrix(isnan(dataMatrix)) = 0;
end

function [columnNames, dataMatrix, phase_numbers] = read_tab_format_without_header(fileName)
    % read_tab_format_without_header Reads a tabular file and separates header, data, and phase count.
    %
    % INPUT:
    %   fileName - Name of the input file (string).
    %
    % OUTPUT:
    %   columnNames - Cell array of column names (header).
    %   dataMatrix - Numeric matrix containing the data.
    %   phase_numbers - Number of mineral phases in the header row.
    %
    % NOTE:
    %   If the file does not contain a valid header or the number of phases is zero,
    %   the function raises an error.

    % Open the file for reading
    fid = fopen(fileName, 'r');
    if fid == -1
        error('Could not open the file: %s', fileName);
    end
    
    % Read the first line to check if it's a header or data
    firstLine = fgetl(fid);
    if isempty(firstLine) || ~ischar(firstLine)
        fclose(fid);
        error('The file %s is empty or has invalid content.', fileName);
    end

    % Check if the first line contains non-numeric characters (indicating a header)
    if any(isstrprop(firstLine, 'alpha'))
        % First line is a header; split it into column names
        columnNames = strsplit(strtrim(firstLine));
        
        % Identify mineral phases (columns beyond the basic ones like node#, P, T, etc.)
        % Assuming phases are all columns starting from the 4th column onward
        phase_numbers = max(0, length(columnNames) - 3);
        
        % If no phases are found, raise an error
        if phase_numbers == 0
            fclose(fid);
            error('The file %s contains no mineral phases. Please ensure the header includes phases.', fileName);
        end
        
        % Read the rest of the file as numeric data
        dataMatrix = [];
        line = fgetl(fid); % Read first data line
        while ischar(line)
            numericRow = sscanf(line, '%f')';
            dataMatrix = [dataMatrix; numericRow];
            line = fgetl(fid);
        end
    else
        % First line is data; no header
        fclose(fid);
        error('The file %s does not contain a valid header. Please ensure the first row contains column names.', fileName);
    end

    % Close the file
    fclose(fid);
end

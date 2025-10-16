function [part1, part2, part3, part4] = extract_mineral_matrices(filename)
    % Read data from file
    fid = fopen(filename, 'r');
    header = fgetl(fid); % Read and skip header line
    data = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f');
    fclose(fid);

    % Extract columns
    depth = data{1};
    pressure_Pa = data{2} * 1e5; % Convert bar to Pa
    temperature = data{3};
    minerals = [data{4}, data{5}, data{6}, data{7}, data{8}, data{9}, data{10}, data{11}, data{12}];
    
    % Mineral indices: O=1, Opx=2, Cpx=3, Gt=4, Wad=5, Ring=6, Pv=7, fp=8, ca-pv=9

    % Part 1: depth < 410km (n*6 matrix)
    idx1 = depth < 410;
    part1 = [depth(idx1), pressure_Pa(idx1), temperature(idx1), ...
             minerals(idx1, 1), ... % O
             minerals(idx1, 2), ... % Opx
             minerals(idx1, 3), ... % Cpx
             minerals(idx1, 4)];    % Gt

    % Part 2: 410km <= depth < 520km (m*4 matrix)
    idx2 = depth >= 410 & depth <= 520;
    part2 = [depth(idx2), pressure_Pa(idx2), temperature(idx2), ...
             minerals(idx2, 5), ... % Wad
             minerals(idx2, 4)];    % Gt

    % Part 3: 520km <= depth <= 660km (k*4 matrix)
    idx3 = depth > 520 & depth <= 660;
    part3 = [depth(idx3), pressure_Pa(idx3), temperature(idx3), ...
             minerals(idx3, 6), ... % Ring
             minerals(idx3, 4)];    % Gt

    % Part 4: depth >= 660km (l*5 matrix)
    idx4 = depth > 660;
    part4 = [depth(idx4), pressure_Pa(idx4), temperature(idx4), ...
             minerals(idx4, 7), ... % Pv
             minerals(idx4, 8), ... % fp
             minerals(idx4, 9)];    % ca-pv

    % Display information
    fprintf('Part 1 (depth < 410km): %d rows x 6 columns\n', size(part1, 1));
    fprintf('Part 2 (410km <= depth < 520km): %d rows x 4 columns\n', size(part2, 1));
    fprintf('Part 3 (520km <= depth <= 660km): %d rows x 4 columns\n', size(part3, 1));
    fprintf('Part 4 (depth >= 660km): %d rows x 5 columns\n', size(part4, 1));
end

function plotModalAbundanceAndCumulativev2(filename)
    % Open the file for reading
    fid = fopen(filename, 'r');
    
    % Skip the first 6 lines of metadata
    for i = 1:6
        fgetl(fid);
    end
    
    % Read the number of rows and columns
    numRows = str2double(fgetl(fid)); % 7th line: number of rows
    numCols = str2double(fgetl(fid)); % 8th line: number of columns
    
    % Read the header line
    headerLine = fgetl(fid);
    headers = strsplit(strtrim(headerLine)); % Extract headers
    
    % Read the data starting from the 9th line
    dataFormat = repmat('%f', 1, numCols); % Format specifier for numeric data
    data = textscan(fid, dataFormat, numRows, 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
    
    % Close the file after reading
    fclose(fid);
    
    % Replace NaN values with zeros in the data
    for i = 1:numel(data)
        data{i}(isnan(data{i})) = 0; % Replace NaN values with zeros
    end
    
    % Extract the x-axis (P(bar)) and y-axis data
    P_bar = data{2}; % Use the second column as x-axis ("P(bar)")
    y_data = cell2mat(data(4:end)); % Use columns 4 onwards as y-axis data, excluding "node#" and "T(K)"
    
    % Calculate cumulative sums for y_data along each row
    cumulative_data = cumsum(y_data, 2); % Compute cumulative sum along columns
    
    % Create a new figure with appropriate size and resolution
    figure('Units', 'inches', 'Position', [0, 0, 12, 12], 'PaperPositionMode', 'auto', 'Color', 'w');
    
    %% Subplot 1: Individual Modal Abundance Curves
    subplot(2, 1, 1); % Top subplot for individual modal abundance
    hold on;
    
    % Use 'lines' colormap for distinct colors for individual modal abundance curves
    distinctColors = lines(size(y_data, 2)); % Generate distinct colors
    for i = 1:size(y_data, 2)
        plot(P_bar / 10000, y_data(:, i), 'LineWidth', 2.0, 'Color', distinctColors(i, :), 'DisplayName', headers{3 + i});
    end
    hold off;
    
    % Set axes properties and title
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'Box', 'on');
    xlabel('P (GPa)', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Modal Abundance (%)', 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'eastoutside', 'FontSize', 12);
    title('Modal Abundance of Minerals vs. P(GPa)', 'FontSize', 18, 'FontWeight', 'bold');
    grid on; % Add grid for better visualization
    
    %% Subplot 2: Cumulative Modal Abundance
    subplot(2, 1, 2); % Bottom subplot for cumulative modal abundance
    hold on;
    
    % Use the manually defined viridis colormap if available
    numColors = size(cumulative_data, 2);
    colors = viridis(numColors); % Use the viridis colormap function

    % Create patch plots for cumulative modal abundance
    for i = 1:numColors
        x_vertices = [P_bar; flipud(P_bar)] / 10000; % X-coordinates
        if i == 1
            y_vertices = [zeros(size(P_bar)); flipud(cumulative_data(:, i))];
        else
            y_vertices = [cumulative_data(:, i - 1); flipud(cumulative_data(:, i))];
        end
        
        % Use colors from the viridis colormap for cumulative patches
        patch(x_vertices, y_vertices, colors(i, :), 'FaceAlpha', 0.6, 'EdgeColor', colors(i, :), 'LineWidth', 1.5);
    end
    
    % Set axes properties and title
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'Box', 'on');
    xlabel('P (GPa)', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Cumulative Modal Abundance (%)', 'FontSize', 16, 'FontWeight', 'bold');
    title('Cumulative Modal Abundance vs. P(GPa)', 'FontSize', 18, 'FontWeight', 'bold');
    xlim([min(P_bar) / 10000, max(P_bar) / 10000]);
    ylim([0, 100]);
    grid on;
    
    % Create a colorbar for the cumulative modal abundance plot
    cbar = colorbar('Ticks', linspace(0, 1, numColors), ...
                    'TickLabels', headers(4:end), ...
                    'Location', 'eastoutside', ...
                    'FontSize', 12); % Remove any duplicate or extra colorbars
                
    % Set colorbar title using the ylabel function (colorbar title cannot be set directly)
    ylabel(cbar, 'Mineral Phases', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Set the colormap axis limits for better visual differentiation
    caxis([0, numColors]); % Adjust color axis limits

    % Remove any additional colorbars that might have been created
    % Remove left or unwanted colorbar if necessary    
    % Save the figure as an EPS file for high-quality printing suitable for Nature publications
    print('plotModalAbundanceAndCumulativeOutput', '-depsc2', '-r300');
    
    % Close the figure to avoid displaying it in the MATLAB window
    close;
end

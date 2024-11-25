    function plotModalAbundanceAndCumulative(filename)
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
    figure('Units', 'inches', 'Position', [0, 0, 12, 6], 'PaperPositionMode', 'auto', 'Color', 'w');
    
    % %% Subplot 1: Individual Modal Abundance Curves
    % subplot(2, 1, 1);
    % hold on;
    % 
    % Use a colormap with faster transitions for better differentiation
    numColors = size(y_data, 2);
    colors = custom_colormap(numColors); % Use custom colormap with rapid transitions
    % for i = 1:numColors
    %     plot(P_bar / 10000, y_data(:, i), 'LineWidth', 2.0, 'Color', colors(i, :), 'DisplayName', headers{3 + i});
    % end
    % hold off;
    % 
    % % Set axes properties and title
    % set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'Box', 'on');
    % xlabel('P (GPa)', 'FontSize', 16, 'FontWeight', 'bold');
    % ylabel('Modal Abundance (%)', 'FontSize', 16, 'FontWeight', 'bold');
    % legend('Location', 'eastoutside', 'FontSize', 12);
    % title('Modal Abundance of Minerals vs. P(GPa)', 'FontSize', 18, 'FontWeight', 'bold');
    % grid on; % Add grid for better visualization
    % 
    % %% Subplot 2: Cumulative Modal Abundance
    % subplot(2, 1, 2);
    % hold on;
    
    % Use patch function to create cumulative modal abundance plot
    for i = 1:size(cumulative_data, 2)
        x_vertices = [P_bar; flipud(P_bar)] / 10000; % X-coordinates
        if i == 1
            y_vertices = [zeros(size(P_bar)); flipud(cumulative_data(:, i))];
        else
            y_vertices = [cumulative_data(:, i - 1); flipud(cumulative_data(:, i))];
        end
        
        % Create filled area plot with transparency and color
        patch(x_vertices, y_vertices, colors(i, :), 'FaceAlpha', 0.7, 'EdgeColor', colors(i, :), 'LineWidth', 1.5);
    end
    
    % Set axes properties and title
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'Box', 'on');
    xlabel('Pressure (GPa)', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Modal Abundance (%)', 'FontSize', 16, 'FontWeight', 'bold');
    % title('Cumulative Modal Abundance vs. P(GPa)', 'FontSize', 18, 'FontWeight', 'bold');
    xlim([min(P_bar) / 10000, max(P_bar) / 10000]);
    ylim([0, 100]);
    % grid on;
    
    % Apply the same custom colormap and add a single colorbar
    colormap(colors); % Apply the custom colormap to the current figure
    cbar = colorbar('Location', 'eastoutside', 'FontSize', 12);
    cbar.Label.String = 'Mineral Phases';
    cbar.Label.FontSize = 16;
    cbar.Label.FontWeight = 'bold';
    
    % Set colorbar ticks and labels to match the headers
    set(cbar, 'Ticks', linspace(0, 1, numColors), 'TickLabels', headers(4:end));
    
    % Save the figure as an EPS file for high-quality printing
    print('plotModalAbundanceAndCumulativeOutput', '-depsc2', '-r300');
    
    % Close the figure to avoid displaying it in the MATLAB window
    close;
end

%% Supporting Function: Custom Colormap with Rapid Transitions
function cmap = custom_colormap(m)
    % Create a custom colormap with 'm' colors and rapid transitions
    base_colors = [
        0.2, 0.4, 1.0; % Blue
        1.0, 0.2, 0.2; % Red
        0.2, 1.0, 0.2; % Green
        1.0, 1.0, 0.2; % Yellow
        1.0, 0.6, 0.2; % Orange
        0.6, 0.2, 1.0; % Purple
        0.2, 1.0, 1.0; % Cyan
        1.0, 0.2, 1.0; % Magenta
        0.8, 0.8, 0.8; % Grey
        0.4, 0.4, 0.4  % Dark Grey
    ];
    
    % Generate the colormap with m colors using interpolation for smooth transitions
    x = linspace(1, size(base_colors, 1), m);
    cmap = interp1(1:size(base_colors, 1), base_colors, x);
end

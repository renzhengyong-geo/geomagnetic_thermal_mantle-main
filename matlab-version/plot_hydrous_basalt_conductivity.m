function plot_hydrous_basalt_conductivity()
    % Reproduce Figure 5(b) from Ni et al. (2011)
    % Electrical conductivity of hydrous basaltic melts at 2 GPa
    
    % Define the VFT equation from the paper
    % logσ = 2.172 - (860.82 - 204.46*w^0.5)/(T - 1146.8)
    
    % Create inverse temperature axis first (main x-axis)
    invT_plot = linspace(5, 7, 100);  % 10000/T range
    T_plot = 10000 ./ invT_plot;      % Convert back to temperature in K
    
    % Water contents (wt%) from the figure legend
    water_contents = [0.02, 1.1, 4.1, 6.3, 0.3];
    labels = {'0.02 wt% H₂O', '1.1 wt% H₂O', '4.1 wt% H₂O', '6.3 wt% H₂O', '0.3 wt% H₂O + 0.5 wt% CO₂'};
    
    % Colors and line styles
    colors = [0 0 1; 0 0.5 0; 1 0 0; 0.5 0 0.5; 0 0 0]; % RGB colors
    line_styles = {'-', '-', '-', '-', '--'};
    line_widths = [2, 2, 2, 2, 2];
    
    % Create figure
    figure('Position', [100, 100, 900, 700]);
    
    % Main axes for conductivity curves
    axes('Position', [0.15, 0.15, 0.7, 0.7]);
    hold on;
    grid on;
    
    % Create cell array to store plot handles for legend
    plot_handles = cell(1, length(water_contents));
       
    % Plot the model prediction curves
    for i = 1:length(water_contents)
        w = water_contents(i);
        conductivity = calculate_conductivity(T_plot, w);
        
        % Plot the predicted curve and store handle
        plot_handles{i} = plot(invT_plot, log10(conductivity), ...
             'Color', colors(i,:), ...
             'LineStyle', line_styles{i}, ...
             'LineWidth', line_widths(i), ...
             'DisplayName', labels{i});
    end
    
    % Set plot properties to match the original figure
    xlabel('10000/T(K)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('log \sigma (S/m)', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Set axis limits to match the original figure
    xlim([5, 7]);
    ylim([0, 2]);
    
    % Set ticks
    set(gca, 'XTick', 5:0.5:7, 'YTick', 0:0.5:2, 'FontSize', 12);
    
    % Add temperature scale on top
    ax1 = gca;
    
    % Create top axis for temperature in °C
    ax2 = axes('Position', ax1.Position, ...
               'XAxisLocation', 'top', ...
               'YAxisLocation', 'right', ...
               'Color', 'none', ...
               'XColor', 'k', 'YColor', 'none');
    
    % Link the x-limits of both axes
    linkaxes([ax1, ax2], 'x');
    
    % Set temperature ticks in °C (correct order: high temp on left, low temp on right)
    temp_ticks_C = [1600, 1500, 1400, 1300, 1200];  % Correct order for inverse temperature
    temp_ticks_K = temp_ticks_C + 273.15;
    invT_ticks_top = 10000 ./ temp_ticks_K;
    
    set(ax2, 'XLim', [5, 7], ...
             'XTick', invT_ticks_top, ...
             'XTickLabel', arrayfun(@(x) sprintf('%d', x), temp_ticks_C, 'UniformOutput', false), ...
             'FontSize', 12);
    xlabel(ax2, 'Temperature (°C)', 'FontSize', 12, 'FontWeight', 'bold');
    
    % Add title
    title('Hydrous Basaltic Melts (2 GPa)', 'FontSize', 14, 'FontWeight', 'bold');
       
    % Create legend using the plot handles
    legend([plot_handles{:}], labels, ...
           'Location', 'northeast', ...
           'FontSize', 11, ...
           'Box', 'on');
    fprintf('Figure reproduced based on Ni et al. (2011) Contrib Mineral Petrol 162:637-650\n');
    fprintf('Curves: Model predictions using VFT equation\n');
end

function sigma = calculate_conductivity(T, w)
    % Calculate electrical conductivity using VFT equation
    % logσ = 2.172 - (860.82 - 204.46*w^0.5)/(T - 1146.8)
    % where:
    % T = temperature in K
    % w = H2O content in wt%
    % sigma = electrical conductivity in S/m
    
    log_sigma = 2.172 - (860.82 - 204.46 * sqrt(w)) ./ (T - 1146.8);
    sigma = 10.^log_sigma;
end
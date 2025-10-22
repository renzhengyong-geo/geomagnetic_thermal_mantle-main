function plot_fO2_relative_to_IW()
    % PLOT_FO2_RELATIVE_TO_IW plots the relative oxygen fugacity
    % (\Delta log(fO2)) of various buffers with respect to the IW buffer.
    %
    % This function uses `compute_fO2` for calculating the absolute fO2 values.

    % It repeats the figure 1 in the following citation
    %  Righter, K., et al. (2023). Oxygen fugacity buffering in high-pressure solid 
    %    media assemblies from IW-6.5 to IW+4.5 and application to the V K-edge oxybarometer.
    %    American Mineralogist, 108(3): 498–513. 
    %    https://doi.org/10.2138/am-2022-8301


    % Buffers to include in the plot
    bufferNames = {'Re_ReO2', 'Ni_NiO', 'Co_CoO', 'W_WO3', 'Fe_FeO', ...
                   'Mo_MoO2', 'Cr_Cr2O3', 'V_V2O3', 'Ta_Ta2O5', ...
                   'Nb_NbO', 'Si_SiO2'};
    
    % Temperature range in Celsius
    T_C = linspace(600, 1700, 500); % Temperature in Celsius

    % Prepare figure
    figure;
    hold on;
    grid on;

    % Line styles and colors for better visibility
    lineStyles = {'-', '--', ':', '-.', '-', '--', ':', '-.', '-', '--', ':'};
    lineWidth = 2.5; % Thicker lines for better visibility
    colorMap = lines(numel(bufferNames)); % Generate distinct colors

    % Compute fO2 for the IW buffer (Fe-FeO)
    IW_fO2 = arrayfun(@(T) compute_fO2_Righter('Fe_FeO', T), T_C);

    % Loop through each buffer and plot relative fO2
    for i = 1:numel(bufferNames)
        bufferName = bufferNames{i};

        % Compute fO2 for the current buffer
        buffer_fO2 = arrayfun(@(T) compute_fO2_Righter(bufferName, T), T_C);

        % Calculate relative fO2 (log10 difference)
        delta_log_fO2 = log10(buffer_fO2) - log10(IW_fO2);

        % Plot the results
        plot(T_C, delta_log_fO2, 'DisplayName', strrep(bufferName, '_', '-'), ...
            'LineStyle', lineStyles{mod(i - 1, numel(lineStyles)) + 1}, ...
            'Color', colorMap(i, :), 'LineWidth', lineWidth);
    end

    % Customize plot
    title('Relative Oxygen Fugacity (\Delta log(fO2)) Buffers', 'FontSize', 14);
    xlabel('Temperature (°C)', 'FontSize', 12);
    ylabel('\Delta log(fO2) relative to IW', 'FontSize', 12);
    legend('Location', 'best');
    set(gca, 'FontSize', 10);

    % Save the figure in high resolution
    exportgraphics(gcf, 'Relative_fO2_Buffers_Plot_Function.png', 'Resolution', 300); % Save as PNG

    % Display the plot
    hold off;
end

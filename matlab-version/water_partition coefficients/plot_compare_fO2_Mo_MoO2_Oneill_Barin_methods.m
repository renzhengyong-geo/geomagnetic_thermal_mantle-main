function plot_compare_fO2_Mo_MoO2_Oneill_Barin_methods()
    % This function compares two methods for calculating oxygen fugacity (fO2)
    % for the Mo-MoO2 buffer: O'Neill (1986) and Barin (1995).
    % 
    % Citations:
    % 1. O'Neill, H. St. C. (1986). Mo-MoO₂ (MOM) oxygen buffer and the free energy 
    %    of formation of MoO₂. American Mineralogist, 71, 1007–1010. 
    %    http://pubs.geoscienceworld.org/msa/ammin/article-pdf/71/7-8/1007/4212801/am71_1007.pdf
    % 
    % 2. Barin, I. (1995). Thermochemical Data of Pure Substances, 3rd ed. 
    %    Wiley-VCH Verlag GmbH.
    % 
    % 3. Righter, K., et al. (2023). Oxygen fugacity buffering in high-pressure solid 
    %    media assemblies from IW-6.5 to IW+4.5 and application to the V K-edge oxybarometer.
    %    American Mineralogist, 108(3): 498–513. 
    %    https://doi.org/10.2138/am-2022-8301

    temperature_C = linspace(0, 2000, 100); % Temperature in Celsius

    % Calculate fO2 using the O'Neill method
    fO2_Oneill = arrayfun(@(T) compute_fO2_Oneill(T), temperature_C);

    % Calculate fO2 using the Barin method (compute_fO2_Righter function)
    fO2_Barin = arrayfun(@(T) compute_fO2_Righter('Mo_MoO2', T), temperature_C);

    % Plot comparison
    figure('Position', [100, 100, 700, 600]);

    % Plot log10(fO2)
    subplot(2, 1, 1);
    plot(temperature_C, log10(fO2_Oneill), 'b--', 'LineWidth', 2, 'DisplayName', 'O''Neill (1986)');
    hold on;
    plot(temperature_C, log10(fO2_Barin), 'r-', 'LineWidth', 2, 'DisplayName', 'Barin (1995)');
    xlabel('Temperature (°C)', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('log_{10}(fO_2)', 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 14);
    title('Comparison of Oxygen Fugacity (fO_2) Methods', 'FontSize', 18, 'FontWeight', 'bold');
    set(gca, 'FontSize', 14, 'LineWidth', 1.5); % Tick marks and axis thickness
    grid on;

    % Plot absolute difference
    subplot(2, 1, 2);
    fO2_diff = abs(fO2_Oneill - fO2_Barin);
    plot(temperature_C, fO2_diff, 'k-', 'LineWidth', 2);
    xlabel('Temperature (°C)', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('|fO_2 Difference|', 'FontSize', 16, 'FontWeight', 'bold');
    title('Absolute Difference Between Methods', 'FontSize', 18, 'FontWeight', 'bold');
    set(gca, 'FontSize', 14, 'LineWidth', 1.5); % Tick marks and axis thickness
    grid on;

    % Save the plot as high-resolution images
    exportgraphics(gcf, 'fO2_Comparison_O_Neill_vs_Barin.png', 'Resolution', 300); % Save as PNG

    hold off;
end

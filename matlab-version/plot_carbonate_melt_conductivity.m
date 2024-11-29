% Temperature range in Celsius
T_Celsius = 950:1:1400;  % From 900°C to 1400°C with a step of 1°C

% Convert temperature to Kelvin
T_Kelvin = T_Celsius + 273.15;

% Calculate the conductivity for each composition
sigma_values = zeros(length(T_Kelvin), 6);  % Pre-allocate matrix for 6 compositions
for i = 1:6
    for j=1:length(T_Kelvin)
        sigma = carbonate_melt_conductivity(T_Kelvin(j));
        sigma_values(j, i) = sigma(i);
    end
end

% Composition labels
composition_labels = {
    '(LiNaK)_2(CO_3)_3', 
    '(NaK)_2(CO_3)_2', 
    '(NaKCa_0.5)_2(CO_3)_3', 
    '(NaKCa)(CO_3)_2', 
    '(KCa_0.5)_2(CO_3)_2', 
    'Mantle Carbonatites'
};

% Plot the results
figure('Units', 'inches', 'Position', [0 0 6 5], 'PaperPositionMode', 'auto');
hold on;
for i = 1:6
    plot(T_Celsius, sigma_values(:, i), 'DisplayName', composition_labels{i}, 'LineWidth', 2);
end

% Formatting the plot
xlabel('Temperature (°C)', 'FontSize', 14);
ylabel('Conductivity (\sigma) [S/m]', 'FontSize', 14);
set(gca, 'YScale', 'log');  % Set y-axis to logarithmic scale
ylim([1e-5, 1e3]);  % Set y-axis limits to match the desired range
xlim([900, 1500]);  % Set y-axis limits to match the desired range

grid on;
legend('show','Location', 'Best');
title('molten carbonate conductivity vs temperature', 'FontSize', 14);
hold off;
box on;

% Save second plot as high-quality image
set(gcf, 'PaperPositionMode', 'auto');
print('molten_carbonate_conductivity_Gaillard_2008', '-dpng', '-r300');

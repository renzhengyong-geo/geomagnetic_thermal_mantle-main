% Plot water fugacity vs. pressure for four different temperatures

% Temperatures in Kelvin to be considered
temperatures = [1273, 1473, 1673, 1873]; % Kelvin

% Pressures to evaluate (in GPa)
pressure_range = linspace(0.1, 15, 50); % Pressure from 0.1 GPa to 15 GPa

% Initialize figure
figure;
hold on;

% Colors for the plot lines
colors = ['b', 'r', 'g', 'm'];

% Loop through each temperature and calculate fugacity for each pressure
for t = 1:length(temperatures)
    temperature_K = temperatures(t);
    fugacity_values = zeros(size(pressure_range));
    
    for p = 1:length(pressure_range)
        pressure_GPa = pressure_range(p);
        % Call the PSfugacity function from water_fugacity_functions
        fugacity_values(p) = water_fugacity_pitzer_sterner_1994('PSfugacity', pressure_GPa, temperature_K - 273.15);
        fugacity_values(p) = water_fugacity_pitzer_sterner(pressure_GPa, temperature_K - 273.15);

    end
    
    % Plot fugacity vs. pressure
    plot(pressure_range, fugacity_values, 'Color', colors(t), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('T = %d K', temperature_K));
end

% Set y-axis to logarithmic scale
set(gca, 'YScale', 'log');
% Set y-axis limits to [10^-5, 10^10]
ylim([1e-5, 1e10]);

% Labels and legend
xlabel('Pressure (GPa)', 'FontSize', 12);
ylabel('Fugacity (GPa)', 'FontSize', 12);
title('Water Fugacity vs Pressure for Different Temperatures', 'FontSize', 14);
legend('show', 'Location', 'NorthWest');
grid on;
box on; 
hold off;

% Save the figure
saveas(gcf, 'water_fugacity_vs_pressure.png'); % Saves the figure as a PNG file

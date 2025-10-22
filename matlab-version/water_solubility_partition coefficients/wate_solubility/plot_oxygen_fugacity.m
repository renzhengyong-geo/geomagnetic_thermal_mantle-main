% Plot log10(fO2) vs. T for a fixed pressure
P_fixed = 5; % GPa
T_range = 1200:50:1800; % Kelvin
log_fO2_T = arrayfun(@(T) oxygen_fugacity(P_fixed, T, 'QFM'), T_range);

figure;
plot(T_range, log_fO2_T, 'LineWidth', 2);
xlabel('Temperature (K)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('log_{10}(fO_2)', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Oxygen Fugacity vs. Temperature (P = %.1f GPa)', P_fixed), ...
    'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);

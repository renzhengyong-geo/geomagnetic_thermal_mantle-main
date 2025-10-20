% CORRECTED testing code for Kohlstedt water solubility model
clear; clc; close all;

% Kohlstedt parameters with proper units
A = 1.1;                    % H/10^6 Si per MPa
n = 1;
delta_V = 10.6e-6;          % m³/mol
R = 8.3145;                 % J/(mol·K)
T_C = 1100;
T_K = T_C + 273.15;

% Pressure range (GPa)
P_GPa = 2:0.5:13;
P_MPa = P_GPa * 1000;       % Convert to MPa

% Realistic fugacity values from Kohlstedt Table 1 (in GPa)
P_data = [2.5, 5.0, 8.0, 9.0, 10.0, 12.0, 13.0];
f_data_GPa = [0.024, 0.070, 2.2, 6.4, 18, 130, 430]; % From Kohlstedt Table 1

% Interpolate fugacity for our pressure range
f_H2O_GPa = interp1(P_data, f_data_GPa, P_GPa, 'pchip', 'extrap');
f_H2O_MPa = f_H2O_GPa * 1000; % Convert to MPa

% Calculate C_OH using Kohlstedt equation
C_OH = A * f_H2O_MPa.^n .* exp(-P_MPa * delta_V / (R * T_K));

% Calculate the ratio (this should DECREASE with pressure!)
C_OH_fH2O_ratio = C_OH ./ f_H2O_MPa;

% Debug output
fprintf('Detailed Calculation:\n');
fprintf('P_GPa\tP_MPa\tf_H2O_MPa\tExpTerm\t\tC_OH\t\tC_OH/f_H2O\n');
fprintf('----------------------------------------------------------------\n');
for i = 1:length(P_GPa)
    exp_term = -P_MPa(i) * delta_V / (R * T_K);
    fprintf('%.1f\t%.0f\t%.1f\t\t%.6f\t%.0f\t\t%.3f\n', ...
            P_GPa(i), P_MPa(i), f_H2O_MPa(i), exp_term, C_OH(i), C_OH_fH2O_ratio(i));
end

% Create plot
figure('Position', [100, 100, 800, 600]);

% Main plot: C_OH/f_H2O vs Pressure (this should show DECREASING trend)
subplot(2,1,1);
plot(P_GPa, C_OH_fH2O_ratio, 'b-', 'LineWidth', 2);
hold on;
plot(P_GPa, C_OH_fH2O_ratio, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
grid on;

xlabel('Pressure (GPa)', 'FontSize', 12);
ylabel('C_{OH} / f_{H_2O} (H/10^6 Si per GPa)', 'FontSize', 12);
title('Kohlstedt et al. (1996) Water Solubility: C_{OH}/f_{H_2O} vs Pressure', 'FontSize', 14);
subtitle('Should show DECREASING trend due to -PΔV/RT term', 'FontSize', 11);

xline(13, '--k', 'α-β Phase Boundary', 'LabelVerticalAlignment', 'bottom');
legend('C_{OH}/f_{H_2O} Ratio', 'Data Points', 'Location', 'northeast');

% Subplot: Log scale to see the exponential decay clearly
subplot(2,1,2);
semilogy(P_GPa, C_OH_fH2O_ratio, 'b-', 'LineWidth', 2);
hold on;
semilogy(P_GPa, C_OH_fH2O_ratio, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
grid on;

xlabel('Pressure (GPa)', 'FontSize', 12);
ylabel('C_{OH} / f_{H_2O} (H/10^6 Si per GPa) - Log Scale', 'FontSize', 12);
title('Exponential Pressure Dependence (Log Scale)', 'FontSize', 12);
xline(13, '--k', 'α-β Phase Boundary', 'LabelVerticalAlignment', 'bottom');

% Calculate theoretical values for comparison
theoretical_ratio = A * exp(-P_MPa * delta_V / (R * T_K));

% Add theoretical curve
subplot(2,1,1);
plot(P_GPa, theoretical_ratio, 'g--', 'LineWidth', 1.5);
legend('C_{OH}/f_{H_2O} Ratio', 'Data Points', 'Theoretical (A*exp(-PΔV/RT))', 'Location', 'northeast');

subplot(2,1,2);
semilogy(P_GPa, theoretical_ratio, 'g--', 'LineWidth', 1.5);
legend('C_{OH}/f_{H_2O} Ratio', 'Data Points', 'Theoretical (A*exp(-PΔV/RT))', 'Location', 'northeast');

% Display results comparison
fprintf('\n=== Comparison with Kohlstedt Table 1 ===\n');
fprintf('P (GPa)\tKohlstedt f_H2O (GPa)\tOur f_H2O (GPa)\tC_OH (H/10^6Si)\n');
fprintf('------------------------------------------------------------\n');

for i = 1:length(P_data)
    idx = find(P_GPa == P_data(i), 1);
    if ~isempty(idx)
        fprintf('%.1f\t%.3f\t\t\t%.3f\t\t%.0f\n', ...
                P_data(i), f_data_GPa(i), f_H2O_GPa(idx), C_OH(idx));
    end
end

% Verify the exponential term
fprintf('\n=== Exponential Term Verification ===\n');
fprintf('At P = 2 GPa: exp(-PΔV/RT) = %.6f\n', exp(-2000 * delta_V / (R * T_K)));
fprintf('At P = 13 GPa: exp(-PΔV/RT) = %.6f\n', exp(-13000 * delta_V / (R * T_K)));
fprintf('Ratio (13GPa/2GPa): %.3f\n', exp(-13000 * delta_V / (R * T_K)) / exp(-2000 * delta_V / (R * T_K)));
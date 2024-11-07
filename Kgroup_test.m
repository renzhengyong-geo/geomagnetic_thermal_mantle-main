clear;

T = 1500;   % Temperature in Kelvin
Cw = logspace(-5, log10(2), 50); % Water content from 10^-5 to 5 (%)
xn = length(Cw);
olivine = zeros(xn, 1);
garnet = zeros(xn, 1);
wadsleyite = zeros(xn, 1);
ringwoodite=zeros(xn,1);
orthopyroxene=zeros(xn,1);
for i = 1:xn
    olivine(i) = Karato_olivine_conductivity(T, Cw(i), 4 * 1e9);
    garnet(i) = Karato_garnet_conductivity(T, Cw(i), 4 * 1e9);
    orthopyroxene(i)=Karato_orthopyroxene_conductivity(T, Cw(i), 4 * 1e9);
    wadsleyite(i) = Karato_wadsleyite_conductivity(T, Cw(i), 15 * 1e9);
    ringwoodite(i)=Karato_ringwoodite_conductivity(T,Cw(i), 15*1e9);
end

% High-quality plot settings
figure('Units', 'inches', 'Position', [0 0 7 6], 'PaperPositionMode', 'auto');

% Plot each mineral's conductivity with specified colors
plot(Cw, olivine, '-', 'LineWidth', 2);
hold on;
plot(Cw, garnet, '-', 'LineWidth', 2);
plot(Cw, orthopyroxene, '-', 'LineWidth', 2);
plot(Cw, wadsleyite, '-', 'LineWidth', 2);
plot(Cw, ringwoodite, '-', 'LineWidth', 2);
hold off;

% Set logarithmic scale for both axes
ylim([1e-5, 10]);
xlim([1e-5, 2]);
set(gca, 'YScale', 'log', 'XScale', 'log', 'FontSize', 12, 'FontName', 'Helvetica');
xlabel('C_w (wt %)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Conductivity (S/m)', 'FontSize', 14, 'FontWeight', 'bold');
title('T = 1500 K', 'FontSize', 16, 'FontWeight', 'bold');

% Adding legend to distinguish between minerals
legend('Olivine', 'Garnet', 'Orthopyroxene','Wadsleyite', 'Ringwoodite' ,'Location', 'best', 'FontSize', 12);

% Customize axes and box
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'TickLength', [0.02 0.02]);
box on;

% Export the plot as a high-resolution image
print(gcf, 'Conductivity_vs_WaterContent', '-dpng', '-r300'); % 300 DPI for publication
close all; % Close any open figures at the start


% This purpose of this script to test the results shown in the citation
% [1] Karato, S. I. (2011). Water distribution across the mantle transition zone
%  and its implications for global material circulation. Earth and Planetary
%  Science Letters, 301(3–4), 413–423. https://doi.org/10.1016/j.epsl.2010.11.038

clear;

T = 1500;   % Temperature in Kelvin
Cw = logspace(-5, log10(2), 50); % Water content from 10^-5 to 5 (%)
xn = length(Cw);
olivine = zeros(xn, 1);
garnet = zeros(xn, 1);
wadsleyite = zeros(xn, 1);
ringwoodite=zeros(xn,1);
orthopyroxene=zeros(xn,1);
clinopyroxene=zeros(xn,1);

for i = 1:xn
    olivine(i) = Karato_olivine_conductivity(T, Cw(i), 4 * 1e9);
    garnet(i) = Karato_garnet_conductivity(T, Cw(i), 4 * 1e9);
    orthopyroxene(i)=Karato_orthopyroxene_conductivity(T, Cw(i), 4 * 1e9);
    clinopyroxene(i)=Karato_clinopyroxene_conductivity(T, Cw(i), 4 * 1e9);
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
plot(Cw, clinopyroxene, '-', 'LineWidth', 2);
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
legend('Olivine', 'Garnet', 'Orthopyroxene', 'Clinopyroxene', 'Wadsleyite', 'Ringwoodite' ,'Location', 'best', 'FontSize', 12);

% Customize axes and box
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'TickLength', [0.02 0.02]);
box on;

% Export the plot as a high-resolution image
print(gcf, 'Conductivity_vs_WaterContent', '-dpng', '-r300'); % 300 DPI for publication
close all; % Close any open figures at the start


% Pyrolite mantle model which is a term used to characterize a model composition of the Earth's mantle. 
% 1) A pyrolitic Upper Mantle is mainly composed of olivine (~60 volume percent (vol%)), clinopyroxene, orthopyroxene, and garnet.
% [7] Pyroxene would gradually dissolved into garnet and form majoritic garnet.[10]

2) A pyrolitic Mantle Transition Zone is mainly composed of 60 vol% olivine-polymorphs (wadsleyite, ringwoodite) and ~40 vol% majoritic garnet. The top and bottom boundary of the Mantle Transition zone are mainly marked by olivine-wadsleyite transition and ringwoodite-perovskite transition, respectively.

3) A pyrolitic Lower Mantle is mainly composed of magnesium perovskite (~80 vol%), ferroperclase (~13 vol%), and calcium perovskite (~7%). In addition, post-perovskite may present at the bottom of the Lower Mantle.
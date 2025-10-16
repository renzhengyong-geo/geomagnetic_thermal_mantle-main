function plot_mineral_abundance(filename)

% Read the entire file
fid = fopen(filename, 'r');
header = fgetl(fid); % Read and skip header line
data = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f');
fclose(fid);

% Extract columns
depth = data{1};
pressure_GPa = data{2} * 1e-4; % Convert bar to GPa
minerals = [data{4:12}]; % 9 minerals

% Get mineral names from header
header_parts = strsplit(header);
mineral_names = header_parts(4:12);

% Create cumulative plot
figure('Position', [100, 100, 600, 600]);
hold on;

% Fixed colors for 9 minerals
colors = [
    0.2  0.2  0.8;   % Dark Blue
    0.8  0.2  0.2;   % Dark Red
    0.2  0.8  0.2;   % Dark Green
    0.8  0.8  0.2;   % Yellow
    0.8  0.2  0.8;   % Magenta
    0.2  0.8  0.8;   % Cyan
    1.0  0.5  0.0;   % Orange
    0.5  0.0  0.5;   % Purple
    0.0  0.5  0.5    % Teal
];

% Create patches for each mineral (from bottom to top)
patch_handles = zeros(1, 9);
for i = 9:-1:1
    if i == 1
        bottom = zeros(size(depth));
        top = minerals(:, 1);
    else
        bottom = sum(minerals(:, 1:i-1), 2);
        top = bottom + minerals(:, i);
    end
    
    x_patch = [depth; flipud(depth)];
    y_patch = [bottom; flipud(top)];
    
    patch_handles(i) = patch(x_patch, y_patch, colors(i,:), ...
        'FaceAlpha', 0.85, 'EdgeColor', 'black', 'LineWidth', 1.0);
end

% Set axes
ylim([0 100]);
xlim([min(depth) max(depth)]);
xlabel('Depth (km)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Modal Abundance (%)', 'FontSize', 14, 'FontWeight', 'bold');

% Create top axis for pressure
ax1 = gca;
ax2 = axes('Position', ax1.Position, ...
           'XAxisLocation', 'top', 'YAxisLocation', 'right', ...
           'Color', 'none', 'XColor', 'r', 'YColor', 'none');

linkaxes([ax1 ax2], 'y');
ax2.XLim = [min(pressure_GPa) max(pressure_GPa)];

% Set ticks using key depths
key_depths = [12.34, 158.48, 287.47, 410, 483.8, 556.7, 660, 809.48];
depth_ticks = arrayfun(@(d) depth(find(depth >= d, 1)), key_depths);
pressure_ticks = arrayfun(@(d) pressure_GPa(find(depth >= d, 1)), key_depths);

% Ensure first and last points
if depth_ticks(1) ~= depth(1)
    depth_ticks = [depth(1), depth_ticks];
    pressure_ticks = [pressure_GPa(1), pressure_ticks];
end
if depth_ticks(end) ~= depth(end)
    depth_ticks = [depth_ticks, depth(end)];
    pressure_ticks = [pressure_ticks, pressure_GPa(end)];
end

set(ax1, 'XTick', depth_ticks);
set(ax2, 'XTick', pressure_ticks);
set(ax2, 'XTickLabel', arrayfun(@(x) sprintf('%.1f', x), pressure_ticks, 'UniformOutput', false));
xlabel(ax2, 'Pressure (GPa)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'r');

% Add legend 
legend(patch_handles, mineral_names, 'Location', 'best', 'FontSize', 11, 'FontWeight', 'bold');

% Improve appearance
grid on;
set(gcf, 'Color', 'white');
set(ax1, 'FontSize', 12, 'FontWeight', 'bold');
set(ax1, 'Box', 'off');
set(ax2, 'Box', 'off');

hold off;

% Display verification
fprintf('First: Depth=%.2f km, Pressure=%.2f GPa\n', depth(1), pressure_GPa(1));
fprintf('Last:  Depth=%.2f km, Pressure=%.2f GPa\n', depth(end), pressure_GPa(end));
fprintf('Row sums: Min=%.1f%%, Max=%.1f%%\n', min(sum(minerals,2)), max(sum(minerals,2)));

% Save the figure as an EPS file for high-quality printing suitable for Nature publications
print('plot_mineral_abundance', '-depsc2', '-r300');

end


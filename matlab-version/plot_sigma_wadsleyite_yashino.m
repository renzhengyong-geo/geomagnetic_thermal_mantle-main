% Define the range for water content and inverse temperature
water_content = linspace(0, 2, 100); % wt% from 0 to 2
inv_temp = linspace(8, 4, 100);      % 10000/T from 8 to 4 (T from 1250K to 2500K, decreasing left to right)

% Preallocate conductivity matrix
conductivity = zeros(length(water_content), length(inv_temp)); % Adjusted dimensions

% Compute conductivity for each combination
for i = 1:length(water_content)
    for j = 1:length(inv_temp)
        T = 10000 / inv_temp(j); % Convert 10000/T back to T (K)
        if T <= 0
            error('Temperature calculated as %f K, which is not greater than zero.', T);
        end
        conductivity(i, j) = Yoshino_wadsleyite_conductivity(T, water_content(i), 0, 0);
    end
end

% Apply logarithmic scale to conductivity (handle zero or negative values)
log_conductivity = log10(max(conductivity, eps)); % Use eps to avoid log(0)

% Create the contour plot with logarithmic scale
figure;
contourf(inv_temp, water_content, log_conductivity, 20, 'LineColor', 'none');
colorbar;
hold on;

% Add contour lines with logarithmic levels matching the plot
C = contour(inv_temp, water_content, log_conductivity, [-2.5 -2 -1.5 -1 -0.5 0 0.5]);
clabel(C, 'FontSize', 10); % Add labels to contours

% Label axes and title
xlabel('10000/T (K)');
ylabel('Water content (wt%)');
title('Conductivity of Wadsleyite (log scale)');

% Add the PEC line (example: linear approximation, adjust based on your data)
PEC_line = 4 + 0.5 * (inv_temp - 4); % Adjust PEC line for decreasing inv_temp from 8 to 4
plot(inv_temp, PEC_line, 'r--', 'LineWidth', 2);

% Set axis limits with increasing order and reverse x-direction
axis([4 8 0 2]); % [xmin xmax ymin ymax] with increasing x
set(gca, 'XDir', 'reverse'); % Reverse x-axis to show 8 to 4 left to right

hold off;
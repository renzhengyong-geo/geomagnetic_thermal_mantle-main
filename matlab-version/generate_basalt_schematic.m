function generate_basalt_schematic()
    % Create figure
    figure('Position', [100, 100, 800, 1000]);
    hold on;
    axis equal;
    
    % Define colors with built-in alpha
    lithosphere_color = [1.0, 0.843, 0.0];  % Gold
    asthenosphere_color = [1.0, 0.7, 0.7];   % Light red
    deep_mantle_color = [1.0, 0.8, 0.6];     % Light orange
    
    % Draw layers
    % Deep mantle
    rectangle('Position', [0, -3, 10, 3], 'FaceColor', deep_mantle_color, 'EdgeColor', 'k');
    % Asthenosphere
    rectangle('Position', [0, 0, 10, 2], 'FaceColor', asthenosphere_color, 'EdgeColor', 'k');
    % Lithosphere
    rectangle('Position', [0, 2, 10, 1], 'FaceColor', lithosphere_color, 'EdgeColor', 'k');
    
    % Draw mid-ocean ridge
    ridge_x = [4.8, 5.0, 5.2];
    ridge_y = [3.0, 2.5, 3.0];
    fill(ridge_x, ridge_y, [0.5, 0.5, 0.5], 'EdgeColor', 'k');
    
    % Draw subduction zone (triangle)
    subduction_x = [1.0, 0.8, 1.2];
    subduction_y = [2.0, 1.5, 1.5];
    fill(subduction_x, subduction_y, [0.65, 0.16, 0.16], 'EdgeColor', 'k');
    
    % Draw hot spots using filled circles
    theta = 0:pi/20:2*pi;
    % Hot spot 1
    x_circle1 = 7 + 0.5*cos(theta);
    y_circle1 = 0 + 0.5*sin(theta);
    fill(x_circle1, y_circle1, 'r', 'EdgeColor', 'r');
    
    % Hot spot 2
    x_circle2 = 8.5 + 0.75*cos(theta);
    y_circle2 = -1 + 0.75*sin(theta);
    fill(x_circle2, y_circle2, [1, 0.6, 0.6], 'EdgeColor', 'r');
    
    % Draw melt pathways with annotation arrows
    % MORB pathway
    annotation('arrow', [0.52, 0.52], [0.35, 0.45], 'Color', 'b', 'LineWidth', 2);
    text(5.2, 1.2, 'MORB', 'Color', 'b', 'FontSize', 11, 'FontWeight', 'bold');
    
    % Subduction zone pathway
    annotation('arrow', [0.12, 0.12], [0.3, 0.4], 'Color', 'g', 'LineWidth', 2);
    text(1.2, 0.8, 'SRB/ARC', 'Color', 'g', 'FontSize', 11, 'FontWeight', 'bold');
    
    % Hot spot pathways
    annotation('arrow', [0.68, 0.68], [0.25, 0.35], 'Color', 'm', 'LineWidth', 2);
    annotation('arrow', [0.82, 0.82], [0.2, 0.3], 'Color', 'm', 'LineWidth', 2);
    text(7.2, -0.3, 'Hot Spot', 'Color', 'm', 'FontSize', 11, 'FontWeight', 'bold');
    
    % OIB pathway
    annotation('arrow', [0.32, 0.32], [0.3, 0.45], 'Color', [1, 0.5, 0], 'LineWidth', 2);
    text(3.2, 1.3, 'OIB', 'Color', [1, 0.5, 0], 'FontSize', 11, 'FontWeight', 'bold');
    
    % Add layer labels
    text(9.3, 2.5, 'Lithosphere', 'FontSize', 12, 'Rotation', 90, 'FontWeight', 'bold');
    text(9.3, 1.0, 'Asthenosphere', 'FontSize', 12, 'Rotation', 90, 'FontWeight', 'bold');
    text(9.3, -1.5, 'Deeper Mantle', 'FontSize', 12, 'Rotation', 90, 'FontWeight', 'bold');
    
    % Add feature labels
    text(4.9, 3.2, 'Mid-Ocean Ridge', 'FontSize', 10, 'HorizontalAlignment', 'center');
    text(1.0, 2.2, 'Subduction Zone', 'FontSize', 10, 'HorizontalAlignment', 'center');
    
    % Set plot properties
    xlim([0, 10]);
    ylim([-3, 4]);
    set(gca, 'XTick', [], 'YTick', []);
    box on;
    
    % Title
    title('Schematic of Basalt Types and Their Mantle Sources', ...
          'FontSize', 14, 'FontWeight', 'bold');
    
    % Add explanation text
    text(0.5, -2.5, 'MORB: Mid-Ocean Ridge Basalt', 'Color', 'b', 'FontSize', 9);
    text(0.5, -2.7, 'SRB/ARC: Subduction-Related Basalt / Arc Basalt', 'Color', 'g', 'FontSize', 9);
    text(0.5, -2.9, 'OIB: Ocean Island Basalt', 'Color', [1, 0.5, 0], 'FontSize', 9);
    
    hold off;
end
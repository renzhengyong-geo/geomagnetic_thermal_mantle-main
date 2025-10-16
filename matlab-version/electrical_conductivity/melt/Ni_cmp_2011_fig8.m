function Ni_cmp_2011_fig8()
    % 复现图8 (from Ni et al. 2011)：地幔电导率随熔体分数的变化
    % 参数设置
    T = 1673;          % 温度 (K)
    D = 0.006;         % 水分配系数
    sigma_mant = 0.01; % 固体地幔电导率 (S/m)
    
    % 源区水含量 (ppm转换为wt%)
    C0_125 = 125 / 10000;  % 125 ppm = 0.0125 wt%
    C0_600 = 600 / 10000;  % 600 ppm = 0.06 wt%
    
    % 熔体分数范围
    F2 = linspace(0.00001, 0.2, 200);  % 熔体分数从0.1%到20%
    
    % 计算两种水含量下的电导率
    sigma_bulk_125 = calculate_bulk_conductivity(F2, T, C0_125, D, sigma_mant);
    sigma_bulk_600 = calculate_bulk_conductivity(F2, T, C0_600, D, sigma_mant);
    
    % 绘图
    figure('Position', [100, 100, 800, 600]);
    
    % 主图 - 对数坐标
    semilogy(F2*100, sigma_bulk_125, 'r-', 'LineWidth', 2, 'DisplayName', '125 ppm H_2O');
    hold on;
    semilogy(F2*100, sigma_bulk_600, 'b-', 'LineWidth', 2, 'DisplayName', '600 ppm H_2O');
    
    xlabel('Melt Fraction (%)');
    ylabel('Bulk Conductivity (S/m)');
    title('Electrical Conductivity vs Melt Fraction (HS+ Model)');
    legend('Location', 'northwest');
    grid on;
    set(gca, 'FontSize', 12);
    
    % 设置Y轴范围以匹配典型的地幔电导率值
    ylim([1e-2, 1]);
    xlim([0, 6]);
    % 添加文本标注参数
    text(0.02, 0.98, sprintf('T = %d K, D = %.3f', T, D), ...
         'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top');
    
    
    % 保存图像
    saveas(gcf, 'mantle_conductivity_vs_melt_fraction.png');
    
    % 显示无水熔体电导率作为参考
    w_dry = 0;
    sigma_melt_dry = calculate_melt_conductivity(w_dry, T);
    fprintf('无水熔体电导率 (w=0) 在 %d K: %.3f S/m\n', T, sigma_melt_dry);
    
    % 显示一些关键点的数值
    fprintf('\n关键点电导率值:\n');
    fprintf('熔体分数\t125 ppm (S/m)\t600 ppm (S/m)\n');
    indices = [1, 25, 50, 75, 100];
    for i = indices
        fprintf('%.1f%%\t\t%.4f\t\t%.4f\n', F2(i)*100, sigma_bulk_125(i), sigma_bulk_600(i));
    end
end

function sigma_bulk = calculate_bulk_conductivity(F2, T, C0, D, sigma_mant)
    % 计算整体电导率
    % 输入:
    %   F2 - 熔体分数数组
    %   T - 温度 (K)
    %   C0 - 源区水含量 (wt%)
    %   D - 水分配系数
    %   sigma_mant - 固体地幔电导率 (S/m)
    % 输出:
    %   sigma_bulk - 整体电导率数组 (S/m)
    
    n = length(F2);
    sigma_bulk = zeros(1, n);
    
    for i = 1:n
        % 1. 计算熔体水含量 (wt%)
        w = C0 / (D + F2(i) * (1 - D));
        
        % 2. 计算熔体电导率 (S/m) - 公式自动处理w=0的情况
        sigma_melt = calculate_melt_conductivity(w, T);
        
        % 3. 使用HS+模型计算整体电导率
        if F2(i) > 0 && F2(i) < 1
            F1 = 1 - F2(i);
            delta_sigma = sigma_melt - sigma_mant;
            
            % HS+ 公式
            denominator = 3 * sigma_melt - F2(i) * delta_sigma;
            if denominator ~= 0
                sigma_bulk(i) = sigma_melt * (1 - (3 * F1 * delta_sigma) / denominator);
            else
                % 避免除零错误
                sigma_bulk(i) = sigma_mant;
            end
        elseif F2(i) >= 1
            sigma_bulk(i) = sigma_melt;
        else
            sigma_bulk(i) = sigma_mant;
        end
    end
end

function sigma_melt = calculate_melt_conductivity(w, T)
    % 计算熔体电导率
    % 输入:
    %   w - 熔体水含量 (wt%)
    %   T - 温度 (K)
    % 输出:
    %   sigma_melt - 熔体电导率 (S/m)
    
    sigma_melt = 10^(2.172 - (860.82 - 204.46 * sqrt(w)) / (T - 1146.8));
end
function plot_Niu_2011_fig9()
    % 复现Niu 2011图9的部分熔融地幔电导率等值线图
    
    % 参数设置
    T_values = [1673, 1573];  % 温度 (K)
    sigma_solid = 0.01;       % 固体地幔电导率 (S/m)
    D = 0.006;                % 水在固-熔体间的分配系数
    
    % 定义网格
    water_content = linspace(3, 6, 100);    % 熔体水含量 3-6 wt%
    melt_fraction = linspace(0.1, 20, 100); % 熔体分数 0.1-20%
    
    % 目标电导率等值线
    target_sigma = [1, 0.3, 0.1, 0.03];  % S/m
    
    % 创建图形
    figure('Position', [100, 100, 800, 400]);
    
    for t_idx = 1:length(T_values)
        T = T_values(t_idx);
        
        % 计算电导率矩阵和对应的地幔源区总水含量
        sigma_bulk = zeros(length(melt_fraction), length(water_content));
        C0_matrix = zeros(length(melt_fraction), length(water_content));
        
        for i = 1:length(melt_fraction)
            F2 = melt_fraction(i) / 100;  % 转换为分数
            F1 = 1 - F2;                  % 固体分数
            
            for j = 1:length(water_content)
                w = water_content(j);     % 熔体水含量 (wt%)
                
                % 计算地幔源区总水含量
                C0 = w * (D + F2 * (1 - D));
                C0_matrix(i,j) = C0;
                
                % 计算熔体电导率
                sigma_melt = calculate_melt_conductivity(w, T);
                
                % HS+ 整体电导率模型
                if F2 > 0 && F2 < 1
                    delta_sigma = sigma_melt - sigma_solid;
                    denominator = 3*sigma_melt - F2*delta_sigma;
                    
                    if denominator ~= 0
                        sigma_bulk(i,j) = sigma_melt * (1 - (3*F1*delta_sigma)/denominator);
                    else
                        sigma_bulk(i,j) = sigma_solid;
                    end
                else
                    sigma_bulk(i,j) = sigma_melt;
                end
            end
        end
        
        % 调试信息
        fprintf('T = %d K: sigma_bulk range = %.4f to %.4f S/m\n', T, min(sigma_bulk(:)), max(sigma_bulk(:)));
        fprintf('T = %d K: C0 range = %.4f to %.4f wt%%\n', T, min(C0_matrix(:)), max(C0_matrix(:)));
        
        % 绘制等值线图
        subplot(1,2,t_idx);
        [C, h] = contour(water_content, melt_fraction, sigma_bulk, target_sigma, 'LineWidth', 2);
        clabel(C, h);
        
        % 设置图形属性
        xlabel('H_2O in Melt (wt%)', 'FontSize', 12);
        ylabel('Melt Fraction (vol%)', 'FontSize', 12);
        title(sprintf('T = %d K', T), 'FontSize', 14);
        
        grid on;
        set(gca, 'YScale', 'log');
        yticks([0.1, 1, 10, 20]);
        yticklabels({'0.1', '1', '10', '20'});
        xlim([3, 6]);
        ylim([0.1, 20]);
    end
    
    % 分析地幔源区总水含量
    analyze_mantle_water_content();
end

function sigma_melt = calculate_melt_conductivity(w, T)
    % 根据Niu 2011论文中的VFT公式计算熔体电导率
    % log σ = 2.172 - (860.82 - 204.46√w) / (T - 1146.8)
    
    numerator = 860.82 - 204.46 * sqrt(w);
    denominator = T - 1146.8;
    log_sigma = 2.172 - (numerator / denominator);
    sigma_melt = 10^log_sigma;
end

function analyze_mantle_water_content()
    % 分析地幔源区总水含量
    fprintf('\n=== 地幔源区总水含量分析 ===\n\n');
    
    D = 0.006;  % 分配系数
    
    % 典型熔体水含量和熔体分数范围
    melt_water_contents = [3, 4, 5, 6];  % wt%
    melt_fractions = [0.1, 1, 5, 10, 20];  % %
    
    fprintf('%-12s %-12s %-15s\n', '熔体水含量', '熔体分数', '地幔源区总水含量');
    fprintf('%-12s %-12s %-15s\n', '(wt%)', '(vol%)', '(wt%)');
    fprintf('%-12s %-12s %-15s\n', '---------', '---------', '---------------');
    
    for w = melt_water_contents
        for F = melt_fractions
            F2 = F / 100;
            % 计算地幔源区总水含量: w = C0 / [D + F2(1 - D)]
            C0 = w * (D + F2 * (1 - D));
            fprintf('%-12.1f %-12.1f %-15.4f\n', w, F, C0);
        end
        fprintf('\n');
    end
    
    % 分析特定电导率条件下的地幔水含量
    fprintf('\n=== 特定电导率条件下的地幔水含量分析 ===\n\n');
    
    T_values = [1573, 1673];
    target_sigma = [1, 0.3, 0.1, 0.03];
    sigma_solid = 0.01;
    
    for T = T_values
        fprintf('温度: %d K\n', T);
        fprintf('%-10s %-12s %-15s %-15s\n', 'σ_target', 'H2O_in_melt', 'Melt_Fraction', 'C0_mantle');
        fprintf('%-10s %-12s %-15s %-15s\n', '(S/m)', '(wt%)', '(vol%)', '(wt%)');
        fprintf('%-10s %-12s %-15s %-15s\n', '--------', '----------', '------------', '----------');
        
        for sigma_target = target_sigma
            % 对于每个熔体水含量，找到达到目标电导率所需的熔体分数
            water_contents = [3, 4, 5, 6];
            
            for w = water_contents
                % 搜索熔体分数
                F2_range = linspace(0.001, 0.2, 1000);
                sigma_calc = zeros(size(F2_range));
                
                for i = 1:length(F2_range)
                    F2 = F2_range(i);
                    sigma_melt = calculate_melt_conductivity(w, T);
                    
                    % HS+ 模型
                    if F2 > 0 && F2 < 1
                        F1 = 1 - F2;
                        delta_sigma = sigma_melt - sigma_solid;
                        denominator = 3*sigma_melt - F2*delta_sigma;
                        if denominator ~= 0
                            sigma_calc(i) = sigma_melt * (1 - (3*F1*delta_sigma)/denominator);
                        else
                            sigma_calc(i) = sigma_solid;
                        end
                    else
                        sigma_calc(i) = sigma_melt;
                    end
                end
                
                % 找到最接近目标电导率的熔体分数
                [~, idx] = min(abs(sigma_calc - sigma_target));
                F2_required = F2_range(idx);
                F_percent = F2_required * 100;
                
                % 计算对应的地幔源区总水含量
                C0_mantle = w * (D + F2_required * (1 - D));
                
                fprintf('%-10.2f %-12.1f %-15.2f %-15.4f\n', ...
                    sigma_target, w, F_percent, C0_mantle);
            end
            fprintf('\n');
        end
        fprintf('----------------------------------------\n\n');
    end
end
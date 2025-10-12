function [F2_HS, F2_parallel, w_HS, w_parallel] = compute_melt_properties(sigma_bulk, C0, T, D, sigma_solid)
    % 使用健壮牛顿迭代法计算给定电导率下的熔体分数和熔体水含量
    
    % 显示输入参数
    fprintf('=== 输入参数 ===\n');
    fprintf('目标整体电导率: %.4f S/m\n', sigma_bulk);
    fprintf('源区总水含量: %.4f wt%%\n', C0);
    fprintf('温度: %.0f K\n', T);
    fprintf('分配系数: %.4f\n', D);
    fprintf('固体电导率: %.4f S/m\n\n', sigma_solid);
    
    % 参数检查
    if T <= 1146.8
        error('温度必须大于1146.8 K');
    end
    if sigma_bulk <= sigma_solid
        error('整体电导率必须大于固体电导率');
    end
    
    % 使用健壮牛顿迭代法求解HS+模型
    fprintf('求解HS+模型...\n');
    [F2_HS, iter_HS, ~] = robust_newton_solver(@(F2) hs_conductivity(F2, C0, T, D, sigma_solid) - sigma_bulk, sigma_bulk, C0, T, 'HS+');
    
    % 使用健壮牛顿迭代法求解平行模型  
    fprintf('求解平行模型...\n');
    [F2_parallel, iter_parallel, ~] = robust_newton_solver(@(F2) parallel_conductivity(F2, C0, T, D, sigma_solid) - sigma_bulk, sigma_bulk, C0, T, 'Parallel');
    
    % 计算对应的熔体水含量
    w_HS = melt_water_content(F2_HS, C0, D);
    w_parallel = melt_water_content(F2_parallel, C0, D);
    
    % 计算熔体电导率
    sigma_melt_HS = melt_conductivity(w_HS, T);
    sigma_melt_parallel = melt_conductivity(w_parallel, T);
    
    % 显示结果
    fprintf('\n=== 计算结果 ===\n');
    fprintf('%-12s %-10s %-12s %-15s %-12s\n', '模型', '熔体分数', '熔体水含量', '熔体电导率', '迭代次数');
    fprintf('%-12s %-10s %-12s %-15s %-12s\n', '------', '--------', '----------', '----------', '--------');
    fprintf('%-12s %-10.2f %-12.2f %-15.2f %-12d\n', ...
        'HS+', F2_HS*100, w_HS, sigma_melt_HS, iter_HS);
    fprintf('%-12s %-10.2f %-12.2f %-15.2f %-12d\n', ...
        'Parallel', F2_parallel*100, w_parallel, sigma_melt_parallel, iter_parallel);
    
    % 绘制结果图形
    plot_results(sigma_bulk, C0, T, D, sigma_solid, F2_HS, F2_parallel, w_HS, w_parallel);
end

function [F2, iter_used, final_residual] = robust_newton_solver(f, sigma_bulk, C0, T, model_name)
    % 健壮牛顿迭代法求解方程 f(F2) = 0
    max_iter = 50;
    tolerance = 1e-10;
    
    % 智能初始猜测
    F2 = get_initial_guess(sigma_bulk, C0, T);
    
    % 跟踪最佳解
    best_F2 = F2;
    best_residual = abs(f(F2));
    
    iter_used = 0;
    
    for iter = 1:max_iter
        iter_used = iter;
        try
            % 计算函数值和数值导数
            [f_val, f_deriv] = compute_derivative(f, F2);
            residual = abs(f_val);
            
            % 更新最佳解
            if residual < best_residual
                best_F2 = F2;
                best_residual = residual;
            end
            
            % 显示迭代信息
            if iter <= 5 || mod(iter, 10) == 0 || iter == max_iter
                fprintf('  %s模型 - 迭代 %2d: F2 = %8.6f, 残差 = %8.2e\n', ...
                    model_name, iter, F2, residual);
            end
            
            % 检查收敛
            if residual < tolerance
                F2 = best_F2;
                final_residual = residual;
                fprintf('  %s模型 - 收敛于迭代 %d\n', model_name, iter);
                return;
            end
            
            % 牛顿迭代步
            if abs(f_deriv) > 1e-15
                delta = f_val / f_deriv;
                
                % 自适应步长控制
                alpha = 1.0;
                max_step = 0.1;
                delta = sign(delta) * min(abs(delta), max_step);
                
                F2_new = F2 - alpha * delta;
                
                % 边界处理
                F2_new = enforce_bounds(F2_new);
                
                % 线搜索确保残差减小
                [F2_new, line_search_success] = line_search(f, F2, F2_new, f_val);
                
                if line_search_success
                    F2 = F2_new;
                else
                    % 线搜索失败，使用备选策略
                    F2 = fallback_strategy(f, F2, sigma_bulk, C0, T, iter);
                end
                
            else
                % 导数为零，使用备选策略
                F2 = fallback_strategy(f, F2, sigma_bulk, C0, T, iter);
            end
            
        catch ME
            % 出错时使用备选策略
            F2 = fallback_strategy(f, F2, sigma_bulk, C0, T, iter);
        end
    end
    
    % 最终收敛检查
    F2 = best_F2;
    final_residual = best_residual;
    if final_residual > 1e-6
        fprintf('  %s模型 - 警告: 未完全收敛，最终残差 = %.2e\n', model_name, final_residual);
    else
        fprintf('  %s模型 - 可接受解，残差 = %.2e\n', model_name, final_residual);
    end
end

function [f_val, f_deriv] = compute_derivative(f, x)
    % 计算函数值和数值导数
    f_val = f(x);
    
    % 自适应步长
    h = max(1e-8, 1e-8 * abs(x));
    
    % 中心差分以获得更精确的导数
    f_plus = f(x + h);
    f_minus = f(x - h);
    f_deriv = (f_plus - f_minus) / (2 * h);
    
    % 如果中心差分有问题，使用前向差分
    if isnan(f_deriv) || isinf(f_deriv)
        f_deriv = (f_plus - f_val) / h;
    end
end

function [x_new, success] = line_search(f, x, x_new, f_val)
    % 线搜索确保残差减小
    max_line_search = 10;
    success = false;
    
    for ls_iter = 1:max_line_search
        try
            f_new = f(x_new);
            
            if abs(f_new) < abs(f_val)
                success = true;
                return;
            end
            
            % 回溯
            x_new = x + 0.5 * (x_new - x);
            
        catch
            % 回溯更多
            x_new = x + 0.25 * (x_new - x);
        end
        
        x_new = enforce_bounds(x_new);
    end
end

function x = fallback_strategy(f, x_current, sigma_bulk, C0, T, iter)
    % 备选策略：当牛顿法遇到问题时使用
    persistent fallback_count
    if isempty(fallback_count)
        fallback_count = 0;
    end
    
    fallback_count = fallback_count + 1;
    
    if fallback_count > 3
        % 多次失败后使用二分法
        x = bisection_search(f, 1e-6, 0.5);
        fallback_count = 0;
    else
        % 尝试新的初始猜测
        x = get_initial_guess(sigma_bulk, C0, T) * (0.8 + 0.4 * rand());
        x = enforce_bounds(x);
    end
end

function x = bisection_search(f, a, b)
    % 二分法搜索
    max_iter = 30;
    tolerance = 1e-8;
    
    fa = f(a);
    fb = f(b);
    
    if fa * fb > 0
        % 区间内无根，使用黄金分割搜索
        x = golden_section_search(f, a, b);
        return;
    end
    
    for iter = 1:max_iter
        x = (a + b) / 2;
        fx = f(x);
        
        if abs(fx) < tolerance || (b - a) < tolerance
            break;
        end
        
        if fa * fx < 0
            b = x;
            fb = fx;
        else
            a = x;
            fa = fx;
        end
    end
end

function x = golden_section_search(f, a, b)
    % 黄金分割搜索找最小值
    golden_ratio = 0.618;
    max_iter = 30;
    tolerance = 1e-8;
    
    for iter = 1:max_iter
        if (b - a) < tolerance
            break;
        end
        
        x1 = b - golden_ratio * (b - a);
        x2 = a + golden_ratio * (b - a);
        
        f1 = abs(f(x1));
        f2 = abs(f(x2));
        
        if f1 < f2
            b = x2;
        else
            a = x1;
        end
    end
    
    x = (a + b) / 2;
end

function x = get_initial_guess(sigma_bulk, C0, T)
    % 智能初始猜测
    if sigma_bulk < 0.05
        x = 0.01;
    elseif sigma_bulk < 0.2
        x = 0.05;
    else
        x = 0.1;
    end
    
    % 根据水含量调整
    if C0 > 0.5
        x = x * 0.8;
    end
    
    % 根据温度调整
    if T > 1600
        x = x * 1.2;
    end
    
    x = enforce_bounds(x);
end

function x = enforce_bounds(x)
    % 强制在合理范围内
    x = max(1e-6, min(0.5, x));
end

function sigma = hs_conductivity(F2, C0, T, D, sigma_solid)
    % 计算HS+模型的整体电导率
    w = melt_water_content(F2, C0, D);
    sigma_melt = melt_conductivity(w, T);
    
    if F2 > 0 && F2 < 1
        F1 = 1 - F2;
        delta_sigma = sigma_melt - sigma_solid;
        denominator = 3*sigma_melt - F2*delta_sigma;
        if abs(denominator) > 1e-12
            sigma = sigma_melt * (1 - (3*F1*delta_sigma)/denominator);
        else
            sigma = sigma_solid;
        end
    else
        sigma = sigma_melt;
    end
end

function sigma = parallel_conductivity(F2, C0, T, D, sigma_solid)
    % 计算平行模型的整体电导率
    w = melt_water_content(F2, C0, D);
    sigma_melt = melt_conductivity(w, T);
    sigma = F2 * sigma_melt + (1 - F2) * sigma_solid;
end

function w = melt_water_content(F2, C0, D)
    % 计算熔体水含量
    w = C0 / (D + F2 * (1 - D));
    w = max(0.1, min(20, w));
end

function sigma_melt = melt_conductivity(w, T)
    % 计算熔体电导率
    log_sigma = 2.172 - (860.82 - 204.46 * sqrt(w)) / (T - 1146.8);
    sigma_melt = 10^log_sigma;
end

function plot_results(sigma_bulk, C0, T, D, sigma_solid, F2_HS, F2_parallel, w_HS, w_parallel)
    % 绘制结果图形
    figure('Position', [100, 100, 1200, 800]);
    
    % 生成数据用于绘图
    F2_range = logspace(-3, -0.3, 200);  % 对数间隔
    sigma_HS = zeros(size(F2_range));
    sigma_parallel = zeros(size(F2_range));
    w_range = zeros(size(F2_range));
    
    for i = 1:length(F2_range)
        F2 = F2_range(i);
        w_range(i) = melt_water_content(F2, C0, D);
        sigma_HS(i) = hs_conductivity(F2, C0, T, D, sigma_solid);
        sigma_parallel(i) = parallel_conductivity(F2, C0, T, D, sigma_solid);
    end
    
    % 子图1: 电导率 vs 熔体分数
    subplot(2,2,1);
    loglog(F2_range*100, sigma_HS, 'b-', 'LineWidth', 2, 'DisplayName', 'HS+ Model');
    hold on;
    loglog(F2_range*100, sigma_parallel, 'r-', 'LineWidth', 2, 'DisplayName', 'Parallel Model');
    yline(sigma_bulk, 'k--', 'LineWidth', 2, 'DisplayName', sprintf('Target: %.3f S/m', sigma_bulk));
    
    % 标记解点
    plot(F2_HS*100, sigma_bulk, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', ...
         'DisplayName', sprintf('HS+: %.2f%%', F2_HS*100));
    plot(F2_parallel*100, sigma_bulk, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
         'DisplayName', sprintf('Parallel: %.2f%%', F2_parallel*100));
    
    xlabel('Melt Fraction (%)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Bulk Conductivity (S/m)', 'FontSize', 12, 'FontWeight', 'bold');
    title('(a) Conductivity vs Melt Fraction', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'northwest');
    grid on;
    set(gca, 'XScale', 'log', 'YScale', 'log');
    
    % 子图2: 熔体水含量 vs 熔体分数
    subplot(2,2,2);
    semilogx(F2_range*100, w_range, 'g-', 'LineWidth', 2, 'DisplayName', 'Water in Melt');
    hold on;
    
    plot(F2_HS*100, w_HS, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', ...
         'DisplayName', sprintf('HS+: %.2f wt%%', w_HS));
    plot(F2_parallel*100, w_parallel, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
         'DisplayName', sprintf('Parallel: %.2f wt%%', w_parallel));
    
    xlabel('Melt Fraction (%)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Water in Melt (wt%)', 'FontSize', 12, 'FontWeight', 'bold');
    title('(b) Water Content vs Melt Fraction', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'northeast');
    grid on;
    set(gca, 'XScale', 'log');
    
    % 子图3: 熔体电导率 vs 熔体分数
    subplot(2,2,3);
    sigma_melt_range = arrayfun(@(w) melt_conductivity(w, T), w_range);
    semilogx(F2_range*100, sigma_melt_range, 'm-', 'LineWidth', 2, 'DisplayName', 'Melt Conductivity');
    
    xlabel('Melt Fraction (%)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Melt Conductivity (S/m)', 'FontSize', 12, 'FontWeight', 'bold');
    title('(c) Melt Conductivity vs Melt Fraction', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    set(gca, 'XScale', 'log', 'YScale', 'log');
    
    % 子图4: 参数摘要
    subplot(2,2,4);
    axis off;
    
    % 创建文本摘要
    text_str = {
        sprintf('输入参数:');
        sprintf('σ_{bulk} = %.3f S/m', sigma_bulk);
        sprintf('C_0 = %.3f wt%%', C0);
        sprintf('T = %d K', T);
        sprintf('D = %.4f', D);
        sprintf('σ_1 = %.3f S/m', sigma_solid);
        '';
        sprintf('HS+ 模型结果:');
        sprintf('F_2 = %.3f%%', F2_HS*100);
        sprintf('w = %.3f wt%%', w_HS);
        '';
        sprintf('平行模型结果:');
        sprintf('F_2 = %.3f%%', F2_parallel*100);
        sprintf('w = %.3f wt%%', w_parallel);
    };
    
    text(0.1, 0.9, text_str, 'FontSize', 12, 'VerticalAlignment', 'top', ...
         'FontName', 'FixedWidth', 'BackgroundColor', [0.95 0.95 0.95]);
    
    title('(d) Parameter Summary', 'FontSize', 14, 'FontWeight', 'bold');
    
    % 添加总标题
    sgtitle(sprintf('Melt Properties Analysis (σ_{bulk} = %.3f S/m, C_0 = %.3f wt%%)', ...
           sigma_bulk, C0), 'FontSize', 16, 'FontWeight', 'bold');
end
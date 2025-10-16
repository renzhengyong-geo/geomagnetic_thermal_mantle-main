function plot_ocean_continental_geotherms()
    % 绘制大洋和大陆地温线
    
    % 深度范围 (km)
    depth = linspace(0, 400, 200);
    
    % 计算地温线
    T_ocean = ocean_geotherm(depth, 50);           % 50 Ma 大洋
    T_cont_stable = continental_geotherm(depth, 50, 15); % 稳定大陆
    T_cont_active = continental_geotherm(depth, 70, 35); % 活动大陆
    
    % 绘图
    figure('Position', [100, 100, 1200, 600]);
    
    % 主图
    subplot(1,2,1);
    plot(T_ocean, depth, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Ocean (50 Ma)');
    hold on;
    plot(T_cont_stable, depth, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Continental (Stable)');
    plot(T_cont_active, depth, 'm-', 'LineWidth', 2.5, 'DisplayName', 'Continental (Active)');
    
    % 绝热参考线
    T_adiabat = adiabat(depth, 1300);
    plot(T_adiabat, depth, 'k--', 'LineWidth', 1, 'DisplayName', 'Adiabat (1300°C)');
    
    set(gca, 'YDir', 'reverse', 'FontSize', 12);
    xlabel('Temperature (°C)', 'FontSize', 14);
    ylabel('Depth (km)', 'FontSize', 14);
    title('Ocean vs Continental Geotherms', 'FontSize', 16);
    legend('Location', 'southeast');
    grid on;
    xlim([0, 1800]);
    ylim([0, 400]);
    
    % 标记重要深度
    yline(410, 'k:', '410 km', 'LabelHorizontalAlignment', 'left');
    yline(660, 'k:', '660 km', 'LabelHorizontalAlignment', 'left');
    
    % 标记LAB
    [~, lab_ocean] = find_lab(T_ocean, depth);
    [~, lab_stable] = find_lab(T_cont_stable, depth);
    [~, lab_active] = find_lab(T_cont_active, depth);
    
    plot(1300, lab_ocean, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    plot(1300, lab_stable, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    plot(1300, lab_active, 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
    
    % 温度梯度图
    subplot(1,2,2);
    
    dTdz_ocean = gradient(T_ocean, depth);
    dTdz_stable = gradient(T_cont_stable, depth);
    dTdz_active = gradient(T_cont_active, depth);
    
    plot(dTdz_ocean, depth, 'b-', 'LineWidth', 2, 'DisplayName', 'Ocean');
    hold on;
    plot(dTdz_stable, depth, 'r-', 'LineWidth', 2, 'DisplayName', 'Stable');
    plot(dTdz_active, depth, 'm-', 'LineWidth', 2, 'DisplayName', 'Active');
    
    set(gca, 'YDir', 'reverse', 'FontSize', 12);
    xlabel('dT/dz (°C/km)', 'FontSize', 14);
    ylabel('Depth (km)', 'FontSize', 14);
    title('Temperature Gradients', 'FontSize', 16);
    legend;
    grid on;
    ylim([0, 400]);
    
    xline(0.5, 'k--', 'Adiabatic', 'LineWidth', 1);
    
    % 输出结果
    fprintf('LAB depths:\n');
    fprintf('  Ocean: %.0f km\n', lab_ocean);
    fprintf('  Stable continental: %.0f km\n', lab_stable);
    fprintf('  Active continental: %.0f km\n', lab_active);
end

function T = ocean_geotherm(z, age)
    % 大洋地温线 - 板块模型 + 绝热延伸
    Tm = 1350;      % 地幔温度 (°C)
    Ts = 0;         % 地表温度 (°C)
    L = 125;        % 板块厚度 (km)
    
    T = zeros(size(z));
    
    for i = 1:length(z)
        if z(i) <= L
            % 板块冷却模型
            T(i) = Ts + (Tm - Ts) * z(i)/L;
        else
            % 绝热延伸
            T_LAB = Ts + (Tm - Ts);
            T(i) = T_LAB + 0.3 * (z(i) - L); % 0.3 K/km
        end
    end
end

function T = continental_geotherm(z, q_s, q_b)
    % 大陆地温线 - 传导 + 对流 + 绝热延伸
    Ts = 10;        % 地表温度 (°C)
    k = 2.5;        % 热导率 (W/m·K)
    A0 = 2e-6;      % 生热率 (W/m³)
    D = 10;         % 生热层厚度 (km)
    Tm = 1300;      % 地幔温度 (°C)
    
    % 转换为标准单位
    q_s = q_s / 1000;   % W/m²
    q_b = q_b / 1000;   % W/m²
    D_m = D * 1000;     % m
    
    % 确定岩石圈厚度
    L = find_L(q_s, Ts, Tm, k, A0, D_m);
    
    T = zeros(size(z));
    
    for i = 1:length(z)
        if z(i) <= L
            % 岩石圈内
            z_m = z(i) * 1000;
            T(i) = Ts + (q_s/k)*z_m + (A0*D_m^2/k)*(1 - exp(-z_m/D_m));
        else
            % 绝热延伸
            z_LAB_m = L * 1000;
            T_LAB = Ts + (q_s/k)*z_LAB_m + (A0*D_m^2/k)*(1 - exp(-L));
            T(i) = T_LAB + 0.3 * (z(i) - L);
        end
    end
end

function L = find_L(q_s, Ts, Tm, k, A0, D)
    % 求解岩石圈厚度
    L_min = 50;
    L_max = 300;
    tolerance = 1;
    
    while (L_max - L_min) > tolerance
        L_mid = (L_min + L_max) / 2;
        L_mid_m = L_mid * 1000;
        
        T_LAB = Ts + (q_s/k)*L_mid_m + (A0*D^2/k)*(1 - exp(-L_mid));
        
        if T_LAB < Tm
            L_min = L_mid;
        else
            L_max = L_mid;
        end
    end
    
    L = (L_min + L_max) / 2;
end

function T = adiabat(z, Tp)
    % 绝热地温线
    T = Tp + 0.3 * z; % 0.3 K/km
end

function [T_lab, z_lab] = find_lab(T, depth)
    % 找到LAB (T = 1300°C)
    [~, idx] = min(abs(T - 1300));
    T_lab = T(idx);
    z_lab = depth(idx);
end
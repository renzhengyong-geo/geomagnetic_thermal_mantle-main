% 示例1: 基本使用
sigma_bulk = 0.1;       % 观测电导率 (S/m)
C0 = 0.0227;            % 源区水含量 227 ppm = 0.0227 wt%
T = 1673;               % 温度 (K)
D = 0.003;              % 分配系数
sigma_solid = 0.01;     % 固体电导率 (S/m)

[F2_HS, F2_parallel, w_HS, w_parallel] = compute_melt_properties(...
    sigma_bulk, C0, T, D, sigma_solid);


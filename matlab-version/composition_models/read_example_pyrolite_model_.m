% 绘图
clear; clc;

% Set output format to long fixed-point (double format)
format short g

% Or use format short g for shorter display
% format short g

filename = 'standard_pyrolite_models.txt';
plot_mineral_abundance(filename);

% 提取矩阵
[part1, part2, part3, part4] = extract_mineral_matrices(filename);

% 查看结果
disp('Part 1 (Olivine, Opx, Cpx, Gt):');
disp(part1);
disp('Part 2 (Wad, Gt):');
disp(part2);
disp('Part 3 (Ring, Gt):');
disp(part3);
disp('Part 4 (Pv, fp, ca-pv):');
disp(part4);
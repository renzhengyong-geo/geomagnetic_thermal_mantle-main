% Example: Two-phase material
f = [0.7, 0.3];           % Volume fractions
sigma = [1.0, 10.0];      % Conductivities

[sigma_upper, sigma_lower] = hashin_shtrikman(f, sigma);

fprintf('Hashin-Shtrikman bounds:\n');
fprintf('Upper bound: %.4f\n', sigma_upper);
fprintf('Lower bound: %.4f\n', sigma_lower);
fprintf('Bounds width: %.4f\n', sigma_upper - sigma_lower);
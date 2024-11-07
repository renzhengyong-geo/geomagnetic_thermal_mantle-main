% Given data
Cw = [0.027, 0.027, 0.080, 0.060, 0.019, 0.019, 0.019, 0.019, 0.013, 0.013, 0.013, 0.013, 0.010, 0.010, 0.010, 0.010, 0.010]; % Water Content (wt%)
T = [960, 1273, 881, 873, 1273, 1073, 973, 873, 1286, 1171, 1068, 963, 1273, 1173, 1073, 973, 873]; % Temperature (K)
s = [0.0025, 0.0386, 0.0010, 0.0009, 0.0273, 0.0032, 0.0014, 0.0004, 0.0133, 0.0055, 0.0026, 0.0007, 0.0108, 0.0060, 0.0036, 0.0017, 0.0004]; % Conductivity (S/m)

% Gas constant R in J/(mol·K)
R = 8.314;

% Define the conductivity model function, where T is in Kelvin (K)
model = @(params, Cw, T) params(1) * Cw.^params(2) .* exp(-params(3)./(R * T));  % Conductivity model

% Initial guess for A, r, and H
initial_guess = [1e-3, 1, 10000];  % Initial estimates for A, r, and H

% Fit the model to the data using nonlinear least squares
options = optimset('Display', 'off');
[params_fit, ~, ~, ~, ~, jacobian] = lsqcurvefit(@(params, Cw) model(params, Cw, T), initial_guess, Cw, s, [], [], options);

% Estimated parameters: A, r, and H
A_fit = params_fit(1);
r_fit = params_fit(2);
H_fit = params_fit(3);

% Covariance matrix and uncertainties
covariance_matrix = inv(jacobian' * jacobian);  % Covariance matrix
uncertainties = sqrt(diag(covariance_matrix));  % Standard deviations (uncertainties)

% Output the results
fprintf('A = %.4e ± %.4e\n', A_fit, uncertainties(1));
fprintf('r = %.4f ± %.4f\n', r_fit, uncertainties(2));
fprintf('H = %.4f ± %.4f\n', H_fit, uncertainties(3));

% Plot the results
Cw_fit = linspace(min(Cw), max(Cw), 100);  % Generating data points for fitting curve
s_fit = model(params_fit, Cw_fit, T(1));   % Evaluate fitted model

figure;
plot(Cw, s, 'o', 'DisplayName', 'Observed Data');  % Plot observed data
hold on;
plot(Cw_fit, s_fit, '-', 'DisplayName', 'Fitted Curve');  % Plot fitted curve
xlabel('Water Content (Cw) [wt%]');
ylabel('Conductivity (s) [S/m]');
legend;
grid on;

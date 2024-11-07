% Given data
Cw = [0.01, 0.019, 0.06, 0.08];   % Water Content (wt%)
T = 873 * ones(size(Cw));          % Fixed Temperature (K)
s = [0.0004, 0.0004, 0.0009, 0.001]; % Conductivity (S/m)

% Gas constant R (J/(mol·K))
R = 8.314;  % J/mol·K

% Logarithmic transformation of the equation
ln_s = log(s);

% Logarithmic transformation of Cw
ln_Cw = log(Cw);

% The exponential term exp(-H/(R*T)) becomes constant since T = 873
exp_term = exp(-1 / (R * 873));   % This is constant for all data points

% Now, modify the equation: ln(s) = ln(A) + r*ln(Cw) - H/(R*T)
% We will fit: ln(s) = ln(A) + r * ln(Cw) - constant_term
X = [ones(length(Cw), 1), ln_Cw']; % Create the design matrix

% Perform least squares fitting to solve for [ln(A), r]
params = X \ (ln_s' - log(exp_term));  % Solve for parameters

% Extract the parameters
ln_A_est = params(1);
r_est = params(2);

% Convert ln(A) back to A
A_est = exp(ln_A_est);

% The constant exponential term was already handled in the fit
H_est = -R * 873 * log(exp_term); % H can be calculated from the constant term

% Display the results
disp(['Estimated A: ', num2str(A_est)]);
disp(['Estimated r: ', num2str(r_est)]);
disp(['Estimated H: ', num2str(H_est)]);

% Plot the data and the fitted model
figure;
loglog(Cw, s, 'bo', 'MarkerFaceColor', 'blue', 'DisplayName', 'Observed Data');
hold on;
s_fit = A_est * Cw.^r_est .* exp(-H_est / (R * 873));  % Predicted conductivity
loglog(Cw, s_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Fitted Model');

% Labels and title
xlabel('Water Content (C_w) [wt %]', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Conductivity (s) [S/m]', 'FontSize', 14, 'FontWeight', 'bold');
title('Fit of Conductivity vs Water Content (T = 873 K)', 'FontSize', 16, 'FontWeight', 'bold');
legend;
grid on;

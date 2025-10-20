function f_H2O = water_fugacity_PitzerSterner(P, T)
% Calculate water fugacity using Pitzer & Sterner (1994) EOS
% Input: P - Pressure in GPa, T - Temperature in C
% Output: f_H2O - Fugacity in GPa

    % Coefficients from Table I
    c_i_j = zeros(10, 6);
    c_i_j(1, :) = [0, 0, 0.24657688e+6, 0.51359951e+2, 0, 0];
    c_i_j(2, :) = [0, 0, 0.58638965e+0, -0.28646939e-2, 0.31375577e-4, 0];  
    c_i_j(3, :) = [0, 0, -0.62783840e+1, 0.14791599e-1, 0.35779579e-3, 0.15432925e-7];
    c_i_j(4, :) = [0, 0, 0, -0.42719875e+0, -0.16325155e-4, 0];
    c_i_j(5, :) = [0, 0, 0.56654978e+4, -0.16580167e+2, 0.76560762e-1, 0];
    c_i_j(6, :) = [0, 0, 0, 0.10917883e+0, 0, 0];
    c_i_j(7, :) = [0.38878656e+13, -0.13494878e+9, 0.30916564e+6, 0.75591105e+1, 0, 0];
    c_i_j(8, :) = [0, 0, -0.65537898e+5, 0.18810675e+3, 0, 0];
    c_i_j(9, :) = [-0.14182435e+14, 0.18165390e+9, -0.19769068e+6, -0.23530318e+2, 0, 0];
    c_i_j(10, :) = [0, 0, 0.92093375e+5, 0.12246777e+3, 0, 0];

    R = 83.14472; % bar·cm3/(mol·K)
    T = T + 273.15; % Convert to K
    P = P * 10000; % Convert GPa to bar
    
    % Calculate molar volume
    rho = 1 / calculate_V(P, T, c_i_j, R);
    
    % Calculate temperature-dependent coefficients c_i(T) (Eq. 4)
    c_i = zeros(1, 10);
    for i = 1:10
        c_i(i) = c_i_j(i,1)*T^(-4) + c_i_j(i,2)*T^(-2) + c_i_j(i,3)*T^(-1) + ...
                 c_i_j(i,4) + c_i_j(i,5)*T + c_i_j(i,6)*T^2;
    end
    
    % Calculate residual Helmholtz energy (Eq. 1)
    % A_res/nRT = c1*rho + [1/(c2 + c3*rho + c4*rho^2 + c5*rho^3 + c6*rho^4) - 1/c2] 
    %            - (c7/c8)[exp(-c8*rho) - 1] - (c9/c10)[exp(-c10*rho) - 1]
    denominator_poly = c_i(2) + c_i(3)*rho + c_i(4)*rho^2 + c_i(5)*rho^3 + c_i(6)*rho^4;
    A_res = c_i(1)*rho + (1/denominator_poly - 1/c_i(2)) ...
            - (c_i(7)/c_i(8))*(exp(-c_i(8)*rho) - 1) ...
            - (c_i(9)/c_i(10))*(exp(-c_i(10)*rho) - 1);
    
    % Calculate compressibility factor Z = P/rhoRT
    Z = P / (rho * R * T);
    
    % Calculate fugacity coefficient phi = exp(A_res + Z - 1 - ln(Z))
    phi = exp(A_res + Z - 1 - log(Z));
    
    % Fugacity f = phi*P
    f_H2O = (phi * P) / 10000; % Convert back to GPa
end

function V = calculate_V(P, T, c_i_j, R)
    % Find molar volume V that satisfies P_eos(V) = P
    P_diff = @(V) P_eos(V, T, c_i_j, R) - P;
    V_guess = max(5.0, R * T / P);
    V = fzero(P_diff, V_guess, optimset('Display', 'off', 'TolX', 1e-12));
end

function P = P_eos(V, T, c_i_j, R)
    % Equation of state (Eq. 2): 
    % P/RT = rho + c1*rho^2 - rho^2[(c3 + 2*c4*rho + 3*c5*rho^2 + 4*c6*rho^3)/(c2 + c3*rho + c4*rho^2 + c5*rho^3 + c6*rho^4)^2] 
    %       + c7*rho^2*exp(-c8*rho) + c9*rho^2*exp(-c10*rho)
    
    rho = 1 / V;
    
    % Calculate temperature-dependent coefficients c_i(T)
    c_i = zeros(1, 10);
    for i = 1:10
        c_i(i) = c_i_j(i,1)*T^(-4) + c_i_j(i,2)*T^(-2) + c_i_j(i,3)*T^(-1) + ...
                 c_i_j(i,4) + c_i_j(i,5)*T + c_i_j(i,6)*T^2;
    end
    
    denominator = c_i(2) + c_i(3)*rho + c_i(4)*rho^2 + c_i(5)*rho^3 + c_i(6)*rho^4;
    numerator = c_i(3) + 2*c_i(4)*rho + 3*c_i(5)*rho^2 + 4*c_i(6)*rho^3;
    
    P_over_RT = rho + c_i(1)*rho^2 - rho^2*(numerator/(denominator^2)) + ...
                c_i(7)*rho^2*exp(-c_i(8)*rho) + c_i(9)*rho^2*exp(-c_i(10)*rho);
    
    P = P_over_RT * R * T;
end
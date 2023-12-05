% Parameters
N = 35;
Q = 1.5;
R = 3;
x0 = 2;
P0 = 8;

% Preallocate arrays
x = zeros(1, N);
y = zeros(1, N);
x(1) = x0 + sqrt(P0) * randn(1);

% Generate state and measurement sequences
for k = 2:N
    x(k) = x(k - 1) + sqrt(Q) * randn(1);
    y(k) = x(k) + sqrt(R) * randn(1);
end

% Incorrect prior at time 0
x0_incorrect = 12; % Change this value to an incorrect prior
P0_incorrect = 8; % Use the same prior covariance as the correct filter

% Initialize arrays with incorrect prior
x_est_incorrect = zeros(1, N);
P_est_incorrect = zeros(1, N);
x_est_incorrect(1) = x0_incorrect;
P_est_incorrect(1) = P0_incorrect;

% Kalman filter loop with incorrect prior
for k = 2:N
    % Time update (prediction)
    x_pred = x_est_incorrect(k - 1);
    P_pred = P_est_incorrect(k - 1) + Q;

    % Measurement update (correction)
    K = P_pred / (P_pred + R); % Kalman gain
    x_est_incorrect(k) = x_pred + K * (y(k) - x_pred);
    P_est_incorrect(k) = (1 - K) * P_pred;
end

% % Plot the results with correct and incorrect priors
% figure;
% % True state sequence
% plot(1:N, x, '-bo', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'True state x_k');
% hold on;
% % Measurement sequence
% plot(1:N, y, 'ro', 'MarkerSize', 6, 'DisplayName', 'Measurement y_k');
% % Estimated state sequence with correct prior
% plot(1:N, x_est, '-go', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Estimated state with correct prior \hat{x}_k|k');
% % Estimated state sequence with incorrect prior
% plot(1:N, x_est_incorrect, '-mo', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Estimated state with incorrect prior \hat{x}_k|k');
% xlabel('Time step k');
% ylabel('Value');
% legend('show');
% grid on;
%% 1D
% Chosen time step k
k_chosen = 7; % Choose a time step that is illustrative

% Gaussian densities for p(xk−1|y1:k−1) and p(xk|y1:k−1)
x_values = linspace(x_est(k_chosen - 1) - 4 * sqrt(P_est(k_chosen - 1)), x_est(k_chosen - 1) + 4 * sqrt(P_est(k_chosen - 1)), 1000);
density_xk1_given_y1_k1 = normpdf(x_values, x_est(k_chosen - 1), sqrt(P_est(k_chosen - 1)));
density_xk_given_y1_k1 = normpdf(x_values, x_est(k_chosen), sqrt(P_est(k_chosen)));

% Gaussian density for p(xk|y1:k)
density_xk_given_y1_k = normpdf(x_values, x_est(k_chosen) + K * (y(k_chosen) - x_est(k_chosen)), sqrt(P_est(k_chosen)));

% Plot the Gaussian densities
figure;
plot(x_values, density_xk1_given_y1_k1, 'b-', 'LineWidth', 2, 'DisplayName', 'p(x_{k-1}|y_{1:k-1})');
hold on;
plot(x_values, density_xk_given_y1_k1, 'r-', 'LineWidth', 2, 'DisplayName', 'p(x_k|y_{1:k-1})');
plot([y(k_chosen), y(k_chosen)], [0, max(density_xk_given_y1_k)], 'g--', 'LineWidth', 2, 'DisplayName', 'y_k');
plot(x_values, density_xk_given_y1_k, 'm-', 'LineWidth', 2, 'DisplayName', 'p(x_k|y_{1:k})');
xlabel('x');
ylabel('Probability Density');
legend('show');
grid on;









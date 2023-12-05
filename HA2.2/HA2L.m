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

% Plot the sequences
figure;
plot(1:N, x, '-bo', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'State sequence x_k');
hold on;
plot(1:N, y, '-ro', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Measurement sequence y_k');
xlabel('Time step k');
ylabel('Value');
legend('show');
grid on;

%% 1b
% Kalman Filter Implementation

% Initialize arrays
x_est = zeros(1, N);
P_est = zeros(1, N);
x_est(1) = x0;
P_est(1) = P0;

% Kalman filter loop
for k = 2:N
    % Time update (prediction)
    x_pred = x_est(k - 1);
    P_pred = P_est(k - 1) + Q;

    % Measurement update (update)
    K = P_pred / (P_pred + R); % Kalman gain
    x_est(k) = x_pred + K * (y(k) - x_pred);
    P_est(k) = (1 - K) * P_pred;
end

% % Plot the results
% figure;
% % True state sequence
% plot(1:N, x, '-bo', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'True state x_k');
% hold on;
% % Measurement sequence
% plot(1:N, y, 'ro', 'MarkerSize', 6, 'DisplayName', 'Measurement y_k');
% % Estimated state sequence
% plot(1:N, x_est, '-go', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Estimated state \hat{x}_k');
% % 3-sigma bounds
% plot(1:N, x_est + 3 * sqrt(P_est), 'g--', 'DisplayName', '3\sigma upper bound');
% plot(1:N, x_est - 3 * sqrt(P_est), 'g--', 'DisplayName', '3\sigma lower bound');
% 
% xlabel('Time step k');
% ylabel('Value');
% legend('show');
% grid on;
%% 1B2
% Time instances to plot
time_instances = [1, 2, 4, 30];

% Plot error density for each specified time instance
figure;
for i = 1:length(time_instances)
    k = time_instances(i);
    err_mean = 0; % Zero-mean
    err_std = sqrt(P_est(k)); % Standard deviation at time k

    % Create an array of error values around the mean
    err_values = linspace(err_mean - 4 * err_std, err_mean + 4 * err_std, 1000);

    % Compute the Gaussian density for each error value
    err_density = normpdf(err_values, err_mean, err_std);

    % Plot the Gaussian density
    subplot(2, 2, i);
    plot(err_values, err_density, 'LineWidth', 2);
    xlabel('Error');
    ylabel('Density');
    title(sprintf('Error Density at k = %d', k));
    grid on;
end





















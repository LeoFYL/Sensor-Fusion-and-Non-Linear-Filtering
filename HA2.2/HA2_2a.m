rng(42); % Set the random seed for reproducibility

n = 35;
Q = 1.5;
R = 3;
x0_mean = 2;
x0_var = 8;

% Generate the state sequence
x = zeros(n, 1);
x(1) = x0_mean + sqrt(x0_var) * randn;
for k = 2:n
    x(k) = x(k-1) + sqrt(Q) * randn;
end

% Generate the measurement sequence
y = x + sqrt(R) * randn(n, 1);

% Plot the state and measurement sequences
figure;
plot(x, 'LineWidth', 2);
hold on;
plot(y, 'o');
hold off;
xlabel('Time step');
ylabel('Value');
title('State and Measurement Sequences');
legend('State sequence (x)', 'Measurement sequence (y)');
grid on;

% Kalman filter initialization
x_prior = x0_mean; % initial state mean
P_prior = x0_var; % initial state covariance

% Pre-allocate arrays for filter output
x_est = zeros(n, 1);
P_est = zeros(n, 1);
x_ub = zeros(n, 1);
x_lb = zeros(n, 1);

% Kalman filter loop
for k = 1:n
    % Predict
    x_pred = x_prior;
    P_pred = P_prior + Q;

    % Update
    K = P_pred / (P_pred + R);
    x_post = x_pred + K * (y(k) - x_pred);
    P_post = (1 - K) * P_pred;

    % Store the results
    x_est(k) = x_post;
    P_est(k) = P_post;
    x_ub(k) = x_post + 3 * sqrt(P_post);
    x_lb(k) = x_post - 3 * sqrt(P_post);

    % Prepare for next iteration
    x_prior = x_post;
    P_prior = P_post;
end

% Plot the estimates with the true states and measurements
figure;
plot(x, 'LineWidth', 2);
hold on;
plot(y, 'o');
plot(x_est, 'r', 'LineWidth', 2);
plot(x_ub, 'r--');
plot(x_lb, 'r--');
hold off;
xlabel('Time step');
ylabel('Value');
title('Kalman Filter Estimates, True States, and Measurements');
legend('True state (x)', 'Measurement (y)', 'Estimate (x\_est)', 'x\_est ± 3σ');
grid on;

% Plot error density for selected time instances
time_instances = [1, 2, 4, 30];
figure;
for i = 1:length(time_instances)
    k = time_instances(i);
    mu = x_est(k) - x(k);
    sigma = sqrt(P_est(k));
    x_range = linspace(mu - 4*sigma, mu + 4*sigma, 100);
    y_density = normpdf(x_range, mu, sigma);
    
    subplot(length(time_instances), 1, i);
    plot(x_range, y_density);
    xlim([-20, 20]);
    xlabel('Error');
    ylabel('Density');
    title(sprintf('Error Density at k = %d', k));
    grid on;
end


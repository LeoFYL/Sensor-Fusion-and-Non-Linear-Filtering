% Parameters
N = 35;
Q = 1.5;
R = 3;
x0 = 2;
P0 = 8;

% Generate a long true state sequence and the corresponding measurement sequence
N_long = 5000; % Length of the sequence
x_long = zeros(1, N_long);
y_long = zeros(1, N_long);
x_long(1) = x0;
y_long(1) = x0 + normrnd(0, sqrt(R));

for k = 2:N_long
    x_long(k) = x_long(k - 1) + normrnd(0, sqrt(Q));
    y_long(k) = x_long(k) + normrnd(0, sqrt(R));
end

% Filter the long measurement sequence using the Kalman filter
x_est_long = zeros(1, N_long);
P_est_long = zeros(1, N_long);
x_est_long(1) = x0;
P_est_long(1) = P0;

for k = 2:N_long
    % Time update (prediction)
    x_pred = x_est_long(k - 1);
    P_pred = P_est_long(k - 1) + Q;

    % Measurement update (update)
    K = P_pred / (P_pred + R); % Kalman gain
    x_est_long(k) = x_pred + K * (y_long(k) - x_pred);
    P_est_long(k) = (1 - K) * P_pred;
end

% Calculate the estimated mean of the sequence
estimated_mean = mean(x_est_long);

% Plot a histogram of the estimation error xk - xˆk|k
estimation_error = x_long - x_est_long;
figure;
histogram(estimation_error, 'Normalization', 'pdf');
hold on;

% Compare the histogram to the pdf N(x; 0, PN|N)
x_values = linspace(-3 * sqrt(P_est_long(N_long)), 3 * sqrt(P_est_long(N_long)), 1000);
density_error = normpdf(x_values, 0, sqrt(P_est_long(N_long)));
plot(x_values, density_error, 'r-', 'LineWidth', 2);
xlabel('Estimation error (x_k - \hat{x}_k|k)');
ylabel('Probability Density');
legend('Histogram', 'N(x; 0, P_N|N)');
grid on;


% Calculate the innovation process vk
innovation = y_long(2:end) - x_est_long(1:end-1);

% Estimate the mean of the innovation process
innovation_mean = mean(innovation);

% Calculate the autocorrelation function of the innovation process
innovation_autocorr = autocorr(innovation);

% Plot the results
figure;

% Plot the autocorrelation function
subplot(2, 1, 1);
stem(innovation_autocorr, 'filled');
xlabel('Lag');
ylabel('Autocorrelation');
title('Autocorrelation Function of the Innovation Process');
grid on;

% Plot the innovation process
subplot(2, 1, 2);
plot(innovation);
xlabel('Time');
ylabel('Innovation');
title(['Innovation Process (Mean = ', num2str(innovation_mean), ')']);
grid on;


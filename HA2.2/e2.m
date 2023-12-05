% Define necessary parameters
x0 = 2;
Q = 1.5;
R = 3;
P0 = 8;
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

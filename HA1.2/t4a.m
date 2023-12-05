% Set the parameters
sigma = 0.5;
p = [0.5, 0.5];
theta_values = [-1, 1];

% Generate a large number of observations of y
num_samples = 10000;
y_samples = zeros(num_samples, 1);
for i = 1:num_samples
    % Sample theta and w
    theta = randsample(theta_values, 1, true, p);
    w = normrnd(0, sigma);
    
    % Compute y
    y = theta + w;
    y_samples(i) = y;
end

% Plot the histogram of y
histogram(y_samples, 'Normalization', 'pdf');
xlabel('y');
ylabel('Probability density');
title('Histogram of y');

% Fit a Gaussian distribution to the data and plot it on top of the histogram
pd = fitdist(y_samples, 'Normal');
x_values = linspace(min(y_samples), max(y_samples), 100);
y_values = pdf(pd, x_values);
hold on;
plot(x_values, y_values, 'r', 'LineWidth', 2);
hold off;
legend('Histogram of y', 'Gaussian fit');

% Define the original distribution of x
mu = 0;
sigma = sqrt(2);
pd = makedist('Normal', mu, sigma);

% Generate a set of random samples from x
n = 10000;
x = random(pd, n, 1);

% Compute the corresponding values of z
z = x.^3;

% Estimate the probability density function of z using kernel density estimation
[f, xi] = ksdensity(z);

% Compute the mean and variance of z
z_mean = mean(z);
z_var = var(z);

% Plot the histogram of the transformed samples and the estimated density of z
figure;
histogram(z, 'Normalization', 'pdf', 'BinWidth', 0.2);
hold on;
plot(xi, f, 'LineWidth', 2);
xlabel('z');
ylabel('Probability density');
title(['PDF of z = x^3, \mu_z = ', num2str(z_mean), ', \sigma_z^2 = ', num2str(z_var)]);
legend('Histogram', 'Estimated density');

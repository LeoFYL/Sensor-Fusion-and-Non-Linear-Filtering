% Parameters
mu_x = 2;
sigma_x = 1;
H = 3;
sigma_r = 0.5;
num_samples = 10000;

% Generate x samples
x = normrnd(mu_x, sigma_x, 1, num_samples);

% Generate r samples
r = normrnd(0, sigma_r, 1, num_samples);

% Compute y samples
y = H*x + r;

% Compute theoretical mean and variance of y
mean_y = H*mu_x;
var_y = H^2*sigma_x^2 + sigma_r^2;

% Compute sample mean and variance of y
sample_mean_y = mean(y);
sample_var_y = var(y);

% Print results
fprintf('Theoretical mean of y: %.2f\n', mean_y)
fprintf('Sample mean of y: %.2f\n', sample_mean_y)
fprintf('Theoretical variance of y: %.2f\n', var_y)
fprintf('Sample variance of y: %.2f\n', sample_var_y)

% Plot histogram of y samples and theoretical normal distribution
histogram(y, 'Normalization', 'pdf')
hold on
x_range = linspace(min(y), max(y), 100);
y_pdf = normpdf(x_range, mean_y, sqrt(var_y));
plot(x_range, y_pdf)
legend('Sample distribution', 'Theoretical distribution')
xlabel('y')
ylabel('pdf')

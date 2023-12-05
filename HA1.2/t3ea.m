


% Simulation parameters
N = 10000;          % number of samples
sigma_r = 1;        % variance of r
a = -1; b = 1;      % interval for x (uniform distribution)

% Generate random variable x with uniform distribution
x = a + (b-a)*rand(N,1);

% Define deterministic functions h(x)
h_linear = @(x) x;
h_nonlinear = @(x) x.^2;

% Calculate y for both linear and non-linear h(x)
y_linear = h_linear(x) + sqrt(sigma_r)*randn(N,1);
y_nonlinear = h_nonlinear(x) + sqrt(sigma_r)*randn(N,1);

% Calculate the true distribution of y
mu_linear = mean(h_linear(x));
mu_nonlinear = mean(h_nonlinear(x));
var_linear = var(h_linear(x))*sigma_r;
var_nonlinear = var(h_nonlinear(x))*sigma_r;
y_min = min([y_linear; y_nonlinear]);
y_max = max([y_linear; y_nonlinear]);
y_grid = linspace(y_min, y_max, 100);
pdf_linear = normpdf(y_grid, mu_linear, sqrt(var_linear));
pdf_nonlinear = normpdf(y_grid, mu_nonlinear, sqrt(var_nonlinear));
pdf_uniform = unifpdf(y_grid, a, b);

% Plot the results
figure;
subplot(2,2,1);
histogram(y_linear, 50, 'Normalization', 'pdf');
hold on;
plot(y_grid, pdf_linear, 'LineWidth', 2);
title('Linear h(x)');
legend('Simulation', 'Theoretical PDF');

subplot(2,2,2);
histogram(y_nonlinear, 50, 'Normalization', 'pdf');
hold on;
plot(y_grid, pdf_nonlinear, 'LineWidth', 2);
title('Non-linear h(x)');
legend('Simulation', 'Theoretical PDF');

subplot(2,2,[3,4]);
plot(y_grid, pdf_linear, 'LineWidth', 2);
hold on;
plot(y_grid, pdf_nonlinear, 'LineWidth', 2);
plot(y_grid, pdf_uniform, 'LineWidth', 2);
title('Distribution of y');
legend('Linear h(x)', 'Non-linear h(x)', 'Uniform x');

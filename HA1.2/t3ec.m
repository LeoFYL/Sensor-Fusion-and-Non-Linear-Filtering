% Simulation parameters
N = 10000;          % number of samples
sigma_r = 1;        % variance of r
a = -1; b = 1;      % interval for x (uniform distribution)

% Generate random variable x with uniform distribution
x = a + (b-a)*rand(N,1);

% Define deterministic function h(x)
H = 2;              % known constant
h = @(x) H*x;

% Calculate y|x and y
y_given_x = h(x) + sqrt(sigma_r)*randn(N,1);
y = h(x) + sqrt(sigma_r)*randn(N,1);

% Calculate the true distribution of y|x and y
mu_given_x = h(x);
var_given_x = sigma_r;
mu_y = H*mean(x);
var_y = H^2*var(x) + sigma_r;

% Generate a grid for plotting
y_min = min([y_given_x; y]);
y_max = max([y_given_x; y]);
y_grid = linspace(y_min, y_max, 100);

% Calculate the true PDFs of y|x and y
pdf_given_x = normpdf(y_grid, mu_given_x, sqrt(var_given_x));
pdf_y = unifpdf(y_grid, mu_y - sqrt(3*var_y), mu_y + sqrt(3*var_y));

% Plot the results
figure;
subplot(1,2,1);
histogram(y_given_x, 50, 'Normalization', 'pdf');
hold on;
plot(y_grid, pdf_given_x, 'LineWidth', 2);
title('Distribution of y given x');
legend('Simulation', 'Theoretical PDF');

subplot(1,2,2);
histogram(y, 50, 'Normalization', 'pdf');
hold on;
plot(y_grid, pdf_y, 'LineWidth', 2);
title('Distribution of y');
legend('Simulation', 'Theoretical PDF');

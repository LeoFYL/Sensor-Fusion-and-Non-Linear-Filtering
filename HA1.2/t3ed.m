% Define parameters
mu_x = 0;
sigma_x = 1;
h = @(x) x.^2; % non-linear function

% Simulate x and r
n_samples = 10000;
x = normrnd(mu_x, sigma_x, n_samples, 1);
r = normrnd(0, sqrt(2), n_samples, 1); % sigma_r^2 = 2

% Calculate y
y = h(x) + r;

% Calculate and plot histograms
figure;
subplot(2,2,1);
histogram(x, 'Normalization', 'pdf');
title('PDF of x');
xlabel('x');
ylabel('PDF');

subplot(2,2,2);
histogram(r, 'Normalization', 'pdf');
title('PDF of r');
xlabel('r');
ylabel('PDF');

subplot(2,2,3);
histogram(y, 'Normalization', 'pdf');
hold on;
xline(mean(y), '--r', 'LineWidth', 1.5);
xline(h(mu_x), '-.b', 'LineWidth', 1.5);
title('PDF of y');
xlabel('y');
ylabel('PDF');
legend('PDF of y', 'Mean of y', 'h(\mu_x)');

% Calculate and plot the true PDFs
x_vals = linspace(mu_x - 5*sigma_x, mu_x + 5*sigma_x, 1000);
y_vals = linspace(h(mu_x) - 5*sqrt(2), h(mu_x) + 5*sqrt(2), 1000);

subplot(2,2,4);
plot(x_vals, normpdf(x_vals, mu_x, sigma_x));
title('True PDF of x');
xlabel('x');
ylabel('PDF');

figure;
plot(y_vals, normpdf(y_vals, h(mu_x), sqrt(2)));
hold on;
histogram(y, 'Normalization', 'pdf');
title('True PDF of y');
xlabel('y');
ylabel('PDF');
legend('True PDF of y', 'PDF of y');

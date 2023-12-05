% Define parameters
mu_x = 0;
sigma_x = sqrt(2);
mu_z = 0;
sigma_z = sqrt(18);

% Define transformation function
transform = @(x) 3*x;

% Generate samples from x
N = 100000;
x_samples = normrnd(mu_x, sigma_x, [N, 1]);

% Transform samples to z
z_samples = transform(x_samples);

% Compute numerical mean and covariance of z
mu_z_num = mean(z_samples);
sigma_z_num = std(z_samples);

% Compute analytical Gaussian pdf of z
z_range = linspace(mu_z - 5*sigma_z, mu_z + 5*sigma_z, 1000);
z_pdf = normpdf(z_range, mu_z, sigma_z);

% Compute numerical pdf of z
[n_counts, z_bins] = hist(z_samples, 100);
z_pdf_num = n_counts / (sum(n_counts) * (z_bins(2) - z_bins(1)));

% Plot results
figure;
hold on;
histogram(z_samples, 100, 'Normalization', 'pdf');
plot(z_range, z_pdf, 'LineWidth', 2);
% plot(z_bins(1:end-1) + (z_bins(2) - z_bins(1))/2, z_pdf_num, 'LineWidth', 2);
xlabel('z');
ylabel('pdf');
legend({'Histogram', 'Analytical Gaussian', 'Numerical Gaussian'});
title('Approximation of Gaussian Transformation');

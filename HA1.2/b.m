% Define probability density function for q
mu_q = [0; 10];
cov_q = [0.3 0; 0 8];
p_q = @(q) mvnpdf(q, mu_q, cov_q);

% Define transformation matrix A
A = [1 0.5; 0 1];

% Calculate mean and covariance of q
mean_q = mu_q;
cov_q = cov_q;

% Calculate mean and covariance of z
mean_z = A * mean_q;
cov_z = A * cov_q * A';

% Plot sigma ellipses for q and z
figure;
hold on;
sigmaEllipse2D(mean_q, cov_q, 1);
sigmaEllipse2D(mean_z, cov_z, 1);
xlabel('q_1');
ylabel('q_2');
title('Mean and Covariance of q and z');
legend('q', 'z');
hold off;

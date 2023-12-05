%% Define the given state densities and sensor locations
mu1 = [125, 125]';
P1 = [100, 0; 0, 25];

mu2 = [-25, 125]';
P2 = [100, 0; 0, 25];

mu3 = [60, 60]';
P3 = [100, 0; 0, 25];

s1 = [0, 100]';
s2 = [100, 0]';

sigma_phi = 0.1 * pi / 180; % Standard deviation of the measurement noise



%% Generate samples for each state density
num_samples = 10000;

samples1 = mvnrnd(mu1, P1, num_samples);
samples2 = mvnrnd(mu2, P2, num_samples);
samples3 = mvnrnd(mu3, P3, num_samples);

%% Compute measurements and add random noise


measurements1 = compute_measurements(samples1, s1, s2, sigma_phi);
measurements2 = compute_measurements(samples2, s1, s2, sigma_phi);
measurements3 = compute_measurements(samples3, s1, s2, sigma_phi);

%% Approximate the mean and covariance of y using the generated samples


[mean1, cov1] = approximate_mean_and_covariance(measurements1);
[mean2, cov2] = approximate_mean_and_covariance(measurements2);
[mean3, cov3] = approximate_mean_and_covariance(measurements3);

% disp('Mean and covariance for state density p1(x):');
% disp('Mean: ');
% disp(mean1);
% disp('Covariance: ');
% disp(cov1);
% 
% disp('Mean and covariance for state density p2(x):');
% disp('Mean: ');
% disp(mean2);
% disp('Covariance: ');
% disp(cov2);
% 
% disp('Mean and covariance for state density p3(x):');
% disp('Mean: ');
% disp(mean3);
% disp('Covariance: ');
% disp(cov3);

%1b
% EKF


[mean1_ekf, cov1_ekf] = ekf_approximation(mu1, P1, s1, s2, sigma_phi);
[mean2_ekf, cov2_ekf] = ekf_approximation(mu2, P2, s1, s2, sigma_phi);
[mean3_ekf, cov3_ekf] = ekf_approximation(mu3, P3, s1, s2, sigma_phi);

% ukf
[mean1_ukf, cov1_ukf] = ukf_approximation(mu1, P1, s1, s2, sigma_phi, 1, 2, 0);
[mean2_ukf, cov2_ukf] = ukf_approximation(mu2, P2, s1, s2, sigma_phi, 1, 2, 0);
[mean3_ukf, cov3_ukf] = ukf_approximation(mu3, P3, s1, s2, sigma_phi, 1, 2, 0);
% CKF
[mean1_ckf, cov1_ckf] = ckf_approximation(mu1, P1, s1, s2, sigma_phi);
[mean2_ckf, cov2_ckf] = ckf_approximation(mu2, P2, s1, s2, sigma_phi);
[mean3_ckf, cov3_ckf] = ckf_approximation(mu3, P3, s1, s2, sigma_phi);


fprintf('\nEKF:\n');
disp('Mean and covariance for state density p1(x):');
disp('Mean: ');
disp(mean1_ekf);
disp('Covariance: ');
disp(cov1_ekf);

disp('Mean and covariance for state density p2(x):');
disp('Mean: ');
disp(mean2_ekf);
disp('Covariance: ');
disp(cov2_ekf);

disp('Mean and covariance for state density p3(x):');
disp('Mean: ');
disp(mean3_ekf);
disp('Covariance: ');
disp(cov3_ekf);

fprintf('\nUKF:\n');
disp('Mean and covariance for state density p1(x):');
disp('Mean: ');
disp(mean1_ukf);
disp('Covariance: ');
disp(cov1_ukf);

disp('Mean and covariance for state density p2(x):');
disp('Mean: ');
disp(mean2_ukf);
disp('Covariance: ');
disp(cov2_ukf);

disp('Mean and covariance for state density p3(x):');
disp('Mean: ');
disp(mean3_ukf);
disp('Covariance: ');
disp(cov3_ukf);

fprintf('\nCKF:\n');
disp('Mean and covariance for state density p1(x):');
disp('Mean: ');
disp(mean1_ckf);
disp('Covariance: ');
disp(cov1_ckf);

disp('Mean and covariance for state density p2(x):');
disp('Mean: ');
disp(mean2_ckf);
disp('Covariance: ');
disp(cov2_ckf);

disp('Mean and covariance for state density p3(x):');
disp('Mean: ');
disp(mean3_ckf);
disp('Covariance: ');
disp(cov3_ckf);










%% Bearing measurement functions
function bearing = bearing_measurement(x, s)
    dx = x - s;
    bearing = atan2(dx(2), dx(1));
end

function dual_bearing = dual_bearing_measurement(x, s1, s2)
    dual_bearing = [bearing_measurement(x, s1); bearing_measurement(x, s2)];
end

function measurements = compute_measurements(samples, s1, s2, sigma_phi)
    measurements = zeros(size(samples, 1), 2);
    noise = sigma_phi * randn(size(samples, 1), 2);

    for i = 1:size(samples, 1)
        measurements(i, :) = dual_bearing_measurement(samples(i, :)', s1, s2)' + noise(i, :);
    end
end

function [mean_y, cov_y] = approximate_mean_and_covariance(measurements)
    mean_y = mean(measurements, 1)';
    cov_y = cov(measurements);
end

%1b

% EKF
function H = compute_jacobian(x, s)
    dx = x - s;
    r = norm(dx);
    H = [(dx(1) / r^2), (dx(2) / r^2)]';
end

function [mean_y, cov_y] = ekf_approximation(mu, P, s1, s2, sigma_phi)
    H1 = compute_jacobian(mu, s1);
    H2 = compute_jacobian(mu, s2);
    H = [H1, H2];
    R = sigma_phi^2 * eye(2);
    
    mean_y = dual_bearing_measurement(mu, s1, s2);
    cov_y = H * P * H' + R;
end

% UKF
function [mean_y, cov_y] = ukf_approximation(mu, P, s1, s2, sigma_phi, alpha, beta, kappa)
    n = length(mu);
    lambda = alpha^2 * (n + kappa) - n;
    wm = [(lambda / (n + lambda)), (1 / (2 * (n + lambda)) * ones(1, 2 * n))];
    wc = [(lambda / (n + lambda) + (1 - alpha^2 + beta)), (1 / (2 * (n + lambda)) * ones(1, 2 * n))];
    R = sigma_phi^2 * eye(2);
    
    % Generate sigma points
    cholP = chol((n + lambda) * P, 'lower');
    sigma_points = [mu, mu + cholP, mu - cholP];
    
    % Transform sigma points
    transformed_points = zeros(2, size(sigma_points, 2));
    for i = 1:size(sigma_points, 2)
        transformed_points(:, i) = dual_bearing_measurement(sigma_points(:, i), s1, s2);
    end
     % Compute mean and covariance
    mean_y = wm * transformed_points';

 cov_y = zeros(2);
for i = 1:size(transformed_points, 2)
    cov_y = cov_y + wc(i) * (transformed_points(:, i) - mean_y) * (transformed_points(:, i) - mean_y)';
end
cov_y = cov_y + R;



end
% CKF
function [mean_y, cov_y] = ckf_approximation(mu, P, s1, s2, sigma_phi)
    n = length(mu);
    R = sigma_phi^2 * eye(2);
    
    % Generate sigma points
    sigma_points = [mu + sqrt(n) * P, mu - sqrt(n) * P];
    
    % Transform sigma points
    transformed_points = zeros(2, size(sigma_points, 2));
    for i = 1:size(sigma_points, 2)
        transformed_points(:, i) = dual_bearing_measurement(sigma_points(:, i), s1, s2);
    end
    
    % Compute mean and covariance
    mean_y = mean(transformed_points, 2);
    cov_y = 1/(2*n) * (transformed_points - mean_y) * (transformed_points - mean_y)' + R;
end




























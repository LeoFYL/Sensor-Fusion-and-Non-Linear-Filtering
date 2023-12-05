% Load sensor measurements
load('SensorMeasurements.mat');

% Design prior (initial state estimate and covariance matrix)
x0 = [0; 0];
P0 = diag([4, 4]);

% Design validation set
n_train = 100;
y_train = Generate_y_seq(n_train);

% Design test set
n_test = 500;
y_test = Generate_y_seq(n_test);

% CV model tuning
sigma_p_cv = 2;
sigma_v_cv = 2;
sigma_p_cv_list = 0.5:0.5:5;
sigma_v_cv_list = 0.5:0.5:5;
rmse_cv = zeros(length(sigma_p_cv_list), length(sigma_v_cv_list));
for i = 1:length(sigma_p_cv_list)
    for j = 1:length(sigma_v_cv_list)
        [state_estimates_cv, rmse_cv(i,j)] = run_kalman_filter(y_train, sigma_p_cv_list(i), sigma_v_cv_list(j), x0, P0, 'cv');
    end
end
[best_i_cv, best_j_cv] = find(rmse_cv == min(rmse_cv(:)));
sigma_p_cv_best = sigma_p_cv_list(best_i_cv);
sigma_v_cv_best = sigma_v_cv_list(best_j_cv);

% CA model tuning
sigma_p_ca = 2;
sigma_v_ca = 2;
sigma_a_ca = 0.5;
sigma_p_ca_list = 0.5:0.5:5;
sigma_v_ca_list = 0.5:0.5:5;
sigma_a_ca_list = 0.1:0.1:1;
rmse_ca = zeros(length(sigma_p_ca_list), length(sigma_v_ca_list), length(sigma_a_ca_list));
for i = 1:length(sigma_p_ca_list)
    for j = 1:length(sigma_v_ca_list)
        for k = 1:length(sigma_a_ca_list)
            [state_estimates_ca, rmse_ca(i,j,k)] = run_kalman_filter_ca(y_train, sigma_p_ca_list(i), sigma_v_ca_list(j), sigma_a_ca_list(k), x0, P0, 'ca');
        end
    end
end
[best_i_ca, best_j_ca, best_k_ca] = find(rmse_ca == min(rmse_ca(:)));
sigma_p_ca_best = sigma_p_ca_list(best_i_ca);
sigma_v_ca_best = sigma_v_ca_list(best_j_ca);
sigma_a_ca_best = sigma_a_ca_list(best_k_ca);

% Plot RMSE surface for CV model
figure;
surf(sigma_p_cv_list, sigma_v_cv_list, rmse_cv');
xlabel('sigma_p');
ylabel('sigma_v');
zlabel('RMSE');
title('RMSE surface for CV model');

% Plot RMSE surface for CA model
figure;
surf(sigma_p_ca_list, sigma_v_ca_list, squeeze(rmse_ca(:, :, best_k_ca)'));
xlabel('sigma_p');
ylabel('sigma_v');
zlabel('RMSE');
title('RMSE surface for CA model');

% Test CV model
[state_estimates_cv, rmse_cv_test] = run_kalman_filter(y_test, sigma_p_cv_best, sigma_v_cv_best, x0, P0, 'cv');
figure;
plot(y_test(1,:), 'b');
hold on;
plot(state_estimates_cv(1,:), 'r');
xlabel('Time step');
ylabel('Position');
title('CV model - Test set');
legend('Measured', 'Estimated');

% Test CA model
[state_estimates_ca, rmse_ca_test] = run_kalman_filter_ca(y_test, sigma_p_ca_best, sigma_v_ca_best, sigma_a_ca_best, x0, P0, 'ca');
figure;
plot(y_test(1,:), 'b');
hold on;
plot(state_estimates_ca(1,:), 'r');
xlabel('Time step');
ylabel('Position');
title('CA model - Test set');
legend('Measured', 'Estimated');

% Print results
fprintf('CV model:\n');
fprintf(' sigma_p: %.2f\n', sigma_p_cv_best);
fprintf(' sigma_v: %.2f\n', sigma_v_cv_best);
fprintf(' RMSE on validation set: %.2f\n', min(rmse_cv(:)));
fprintf(' RMSE on test set: %.2f\n', rmse_cv_test);
fprintf('CA model:\n');
fprintf(' sigma_p: %.2f\n', sigma_p_ca_best);
fprintf(' sigma_v: %.2f\n', sigma_v_ca_best);
fprintf(' sigma_a: %.2f\n', sigma_a_ca_best);
fprintf(' RMSE on validation set: %.2f\n', min(rmse_ca(:)));
fprintf(' RMSE on test set: %.2f\n', rmse_ca_test);

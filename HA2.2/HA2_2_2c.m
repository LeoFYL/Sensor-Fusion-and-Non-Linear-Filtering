% Load sensor data using Generate_y_seq() function
y_seq = Generate_y_seq();

% Split the sensor data into a training set (80%) and a validation set (20%)
train_frac = 0.8;
n_train = round(size(y_seq, 2) * train_frac);
y_train = y_seq(:, 1:n_train);
y_val = y_seq(:, n_train+1:end);

% Initialize the search space for position measurement noise variance (R_p)
R_p_values = linspace(1, 4, 10);

% Initialize the search space for process noise covariance matrix (Q) elements
sigma_values = linspace(0.1, 10, 10);

% Initialize variables to store the best parameters and the minimum RMSE
best_R_p = 0;
best_sigma = [0, 0];
min_rmse = Inf;

% Loop through all combinations of R_p and sigma values
for R_p = R_p_values
    for sigma_p = sigma_values
        for sigma_v = sigma_values
            % Run the Kalman filter for the CV model with the current parameters
            [state_estimates, rmse] = run_kalman_filter(y_train, R_p, sigma_p, sigma_v, 'cv');
            
            % Update the best parameters if the current RMSE is smaller
            if rmse < min_rmse
                min_rmse = rmse;
                best_R_p = R_p;
                best_sigma = [sigma_p, sigma_v];
            end
        end
    end
end
% Run the Kalman filter for the CV model with the best parameters on the validation set
[state_estimates_cv, rmse_cv] = run_kalman_filter(y_val, best_R_p, best_sigma(1), best_sigma(2), 'cv');

% Plot the state estimates for the CV model
figure;
subplot(2, 1, 1);
plot(state_estimates_cv(1, :));
title('Estimated Position (CV Model)');
xlabel('Time step');
ylabel('Position (m)');

subplot(2, 1, 2);
plot(state_estimates_cv(2, :));
title('Estimated Speed (CV Model)');
xlabel('Time step');
ylabel('Speed (m/s)');
% Define run_kalman_filter function
function [state_estimates, rmse] = run_kalman_filter(y_data, R_p, sigma_p, sigma_v, model_type)
    % Set up the Kalman filter according to the given model type (CV or CA)
    if strcmp(model_type, 'cv')
        % Time step
        dt = 0.1;

        % State transition matrix
        A = [1, dt; 0, 1];

        % Process noise covariance matrix
        Q = diag([sigma_p, sigma_v]);

        % Measurement matrix
        H = [1, 0];
    end
    
    % Initialization
    x = [0; 0]; % Initial state [position; speed]
    P = eye(2); % Initial state covariance

    % Initialize variables to store the state estimates and RMSE
    state_estimates = zeros(2, size(y_data, 2));
    rmse = 0;

    % Loop through the sensor data
    for k = 1:size(y_data, 2)
        % Prediction step
        x = A * x;
        P = A * P * A' + Q;

        % If there is a position measurement available, perform the update step
        if ~isnan(y_data(1, k))
            % Measurement update step
            K = P * H' / (H * P * H' + R_p);
            x = x + K * (y_data(1, k) - H * x);
            P = (eye(2) - K * H) * P;
        end

        % Store the state estimate for the current time step
        state_estimates(:, k) = x;

        % Calculate the RMSE if the true position is available (for tuning purposes)
        if ~isnan(y_data(2, k))

            rmse = rmse + (y_data(2, k) - x(1))^2;
        end
    end

    % Calculate the final RMSE
    rmse = sqrt(rmse / size(y_data, 2));
end
















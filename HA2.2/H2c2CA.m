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
best_sigma = [0, 0, 0];
min_rmse = Inf;

% Loop through all combinations of R_p and sigma values
for R_p = R_p_values
    for sigma_p = sigma_values
        for sigma_v = sigma_values
            for sigma_a = sigma_values
                % Run the Kalman filter for the CA model with the current parameters
                [state_estimates, rmse] = run_kalman_filter(y_train, R_p, sigma_p, sigma_v, sigma_a, 'ca');
                
                % Update the best parameters if the current RMSE is smaller
                if rmse < min_rmse
                    min_rmse = rmse;
                    best_R_p = R_p;
                    best_sigma = [sigma_p, sigma_v, sigma_a];
                end
            end
        end
    end
end
% Run the Kalman filter for the CA model with the best parameters on the validation set
[state_estimates_ca, rmse_ca] = run_kalman_filter(y_val, best_R_p, best_sigma(1), best_sigma(2), best_sigma(3), 'ca');

% Plot the state estimates for the CA model
figure;
subplot(3, 1, 1);
plot(state_estimates_ca(1, :));
title('Estimated Position (CA Model)');
xlabel('Time step');
ylabel('Position (m)');

subplot(3, 1, 2);
plot(state_estimates_ca(2, :));
title('Estimated Speed (CA Model)');
xlabel('Time step');
ylabel('Speed (m/s)');

subplot(3, 1, 3);
plot(state_estimates_ca(3, :));
title('Estimated Acceleration (CA Model)');
xlabel('Time step');
ylabel('Acceleration (m/s^2)');

% Define run_kalman_filter function
function [state_estimates, rmse] = run_kalman_filter(y_data, R_p, sigma_p, sigma_v, sigma_a, model_type)
    % Set up the Kalman filter according to the given model type (CV or CA)
    if strcmp(model_type, 'ca')
        % Time step
        dt = 0.1;

        % State transition matrix
        A = [1, dt, dt^2/2; 0, 1, dt; 0, 0, 1];

        % Process noise covariance matrix
        Q = diag([sigma_p, sigma_v, sigma_a]);

        % Measurement matrix
        H = [1, 0, 0];
    end
    
    % Initialization
    x = [0; 0; 0]; % Initial state

% (m) [position; speed; acceleration]
P = eye(3); % Initial state covariance

% Initialize variables to store the state estimates and RMSE
state_estimates = zeros(3, size(y_data, 2));
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
        P = (eye(3) - K * H) * P;
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































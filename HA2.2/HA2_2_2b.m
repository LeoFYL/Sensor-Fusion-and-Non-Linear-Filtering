% Load sensor data using Generate_y_seq() function
% Assuming the output data is stored in a matrix y_seq
% with position measurements in the first row and speed measurements in the second row
y_seq = Generate_y_seq();

% Kalman filter initialization
% (Adjust these initial values according to your specific problem)
x = [0; 0]; % Initial state estimate [position; speed]
P = eye(2); % Initial state covariance matrix
Q = diag([1, 1]); % Process noise covariance matrix
R_p = 1; % Measurement noise covariance for position
R_v = 1; % Measurement noise covariance for speed

% Time step
dt = 0.1;

% State transition matrix
A = [1, dt; 0, 1];

% Measurement matrices for position and speed
H_p = [1, 0];
H_v = [0, 1];

% Kalman filter loop
n = size(y_seq, 2);
state_estimates = zeros(2, n);

for i = 1:n
    % Prediction step
    x = A * x;
    P = A * P * A' + Q;

    % Check if position measurement is available
    if ~isnan(y_seq(1, i))
        % Update step with both position and speed measurements
        z = y_seq(:, i); % Measurement vector
        H = [H_p; H_v]; % Combined measurement matrix
        R = diag([R_p, R_v]); % Combined measurement noise covariance matrix
        
        % Kalman gain
        K = P * H' / (H * P * H' + R);
        
        % Update state and covariance
        x = x + K * (z - H * x);
        P = (eye(2) - K * H) * P;
    else
        % Update step with only speed measurement
        z = y_seq(2, i); % Measurement vector
        H = H_v; % Measurement matrix
        R = R_v; % Measurement noise covariance matrix
        
        % Kalman gain
        K = P * H' / (H * P * H' + R);
        
        % Update state and covariance
        x = x + K * (z - H * x);
        P = (eye(2) - K * H) * P;
    end
    
    % Store state estimates
    state_estimates(:, i) = x;
end

% Plot the state estimates
figure;
subplot(2, 1, 1);
plot(state_estimates(1, :));
title('Estimated Position');
xlabel('Time step');
ylabel('Position (m)');

subplot(2, 1, 2);
plot(state_estimates(2, :));
title('Estimated Speed');
xlabel('Time step');
ylabel('Speed (m/s)');

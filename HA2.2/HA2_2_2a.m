

load('SensorMeasurements.mat');

% Calculate the mean of measured speed values for each trial
stationary_mean = mean(CalibrationSequenceVelocity_v0);
ten_mps_mean = mean(CalibrationSequenceVelocity_v10);
twenty_mps_mean = mean(CalibrationSequenceVelocity_v20);

% Known true speeds
true_speeds = [0, 10, 20];

% Calculate scaling constant C for each trial and average them
C_estimates = [stationary_mean, ten_mps_mean, twenty_mps_mean] ./ true_speeds;
C = nanmean(C_estimates);  % Use nanmean to ignore NaN values (when dividing by 0)

% Calculate the noise r_k^v for each trial
stationary_noise = (CalibrationSequenceVelocity_v0 / C) - true_speeds(1);
ten_mps_noise = (CalibrationSequenceVelocity_v10 / C) - true_speeds(2);
twenty_mps_noise = (CalibrationSequenceVelocity_v20 / C) - true_speeds(3);

% Calculate the variance of the noise for each trial and average them
Var_rkv = mean([var(stationary_noise), var(ten_mps_noise), var(twenty_mps_noise)]);

% Display the results
fprintf('Estimated scaling constant C: %f\n', C);
fprintf('Estimated variance of the velocity sensor noise Var[r_k^v]: %f\n', Var_rkv);

% Plot actual speed vs. measured speed before calibration
figure;
subplot(2,1,1);
plot(CalibrationSequenceVelocity_v0, 'o');
hold on;
plot(CalibrationSequenceVelocity_v10, 'x');
plot(CalibrationSequenceVelocity_v20, 's');
xlabel('Sample index');
ylabel('Measured speed (m/s)');
title('Actual speed vs. measured speed (Before Calibration)');
legend('0 m/s', '10 m/s', '20 m/s');
hold off;

% Plot actual speed vs. measured speed after calibration
subplot(2,1,2);
plot(CalibrationSequenceVelocity_v0 / C, 'o');
hold on;
plot(CalibrationSequenceVelocity_v10 / C, 'x');
plot(CalibrationSequenceVelocity_v20 / C, 's');
xlabel('Sample index');
ylabel('Measured speed (m/s)');
title('Actual speed vs. measured speed (After Calibration)');
legend('0 m/s', '10 m/s', '20 m/s');
hold off;

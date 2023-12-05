function Y = genLinearMeasurementSequence(X, H, R)
%GENLINEARMEASUREMENTSEQUENCE generates a sequence of observations of the state 
% sequence X using a linear measurement model. Measurement noise is assumed to be 
% zero mean and Gaussian.
%
%Input:
%   X           [n x N+1] State vector sequence. The k:th state vector is X(:,k+1)
%   H           [m x n] Measurement matrix
%   R           [m x m] Measurement noise covariance
%
%Output:
%   Y           [m x N] Measurement sequence
%
 % Get dimensions
    [n, N_plus_1] = size(X);
    [m, n_H] = size(H);
    
    % Validate input dimensions
    if n ~= n_H
        error('H matrix dimensions do not match with state sequence');
    end
    if m ~= size(R, 1) || m ~= size(R, 2)
        error('R matrix dimensions do not match with measurement matrix H');
    end

    % Preallocate memory for the measurement sequence
    Y = zeros(m, N_plus_1 - 1);

    % Generate measurement sequence using the linear measurement model
    for k = 1:(N_plus_1 - 1)
        Y(:, k) = H * X(:, k+1) + mvnrnd(zeros(m, 1)', R)';
    end
% your code here
end
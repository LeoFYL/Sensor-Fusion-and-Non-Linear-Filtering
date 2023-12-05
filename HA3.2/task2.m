%% 2.2 a
clc
rng(970427);

Q = 1.5;
R = 2.5;
x0 = 2;
P0 = 6;

N = 100; %Number of particles

A = 1;
f = @(x) A*x; % Process Model
C = 1; 
h = @(x) C*x; % Measurement Model
T = 0.1; % Sampling time
K = 2/T; % 2s long data

X = genLinearStateSequence(x0, P0, A, Q, K);
Y = genLinearMeasurementSequence(X, C, R);

[X_kal, P_kal] = kalmanFilter(Y, x0, P0, A, Q, C, R);
resampling = false;
kernal = 1;
plotFunc_handle  = @(k, Xk, Xkmin1, Wk, j) ...
                    (plotPostPdf(k, Xk, Wk, X_kal, P_kal, resampling, kernal));
[xfp, Pfp, Xp, Wp] = pfFilter(x0, P0, Y, f, Q, h, R, N, resampling, plotFunc_handle); %Plot func not in arguments
resampling = true;
plotFunc_handle  = @(k, Xk, Xkmin1, Wk, j) ...
                    plotPostPdf(k, Xk, Wk, X_kal, P_kal, resampling, kernal);
[xfp_resamp, Pfp_resamp, Xp_resamp, Wp_resamp] = pfFilter(x0, P0, Y, f, Q, h, R, N, resampling, plotFunc_handle); %Plot func not in arguments




figure(1); clf; hold on; grid on;
plot((0:K).*T ,X, 'b', 'LineWidth', 2)              % True state
plot((1:K).*T,Y, '*m')                              % Measurements
plot((0:K).*T, [x0 X_kal], '-r', 'LineWidth', 2)    % Kalman filter estimate
var = P_kal(:,:,:);                                 % Kalman 3 sigma    
plot((0:K).*T, [x0 X_kal] + 3 ...
    *sqrt([P0 var(:)']),'--', 'color', '#f07382');
plot((0:K).*T, [x0 X_kal] - 3*sqrt([P0 var(:)']), '--', 'color', '#f07382','HandleVisibility','off');

plot((0:K).*T, [x0 xfp], '-g', 'LineWidth', 2)     % Particle filter estimate without resampling   
var = Pfp(:,:,:);                                   % Particle 3 sigma without resampling
plot((0:K).*T, [x0 xfp] + 3*sqrt([P0 var(:)']),'--', 'color', '#74c971');
plot((0:K).*T, [x0 xfp] - 3*sqrt([P0 var(:)']), '--', 'color', '#74c971','HandleVisibility','off');


plot((0:K).*T, [x0 xfp_resamp], '-m', 'LineWidth', 2)     % Particle filter estimate with resampling  
var = Pfp_resamp(:,:,:);                                   % Parficle 3 sigma with resampling
plot((0:K).*T, [x0 xfp_resamp] + 3*sqrt([P0 var(:)']),'--', 'color', '#d687d2');
plot((0:K).*T, [x0 xfp_resamp] - 3*sqrt([P0 var(:)']), '--', 'color', '#d687d2', 'HandleVisibility','off');

title('Kalman vs Particle filter')
xlabel('Time [s]');
ylabel('State value');
legend('True state ', 'Measurements', 'Kalman','$\pm 3\sigma$ Kalman', ...
    'PF w/o resampling', '$\pm 3\sigma$ PF w/o resampling',  ...
    'PF w/ resampling', '$\pm 3\sigma$ PF w/ resampling', 'Interpreter','latex')



mse.kalman = sum((X(2:end)-X_kal).^2);
mse.pf = sum((X(2:end)-xfp).^2);







function [xs, Ps, xf, Pf, xp, Pp] = ...
    nonLinRTSsmoother(Y, x_0, P_0, f, T, Q, S, h, R, sigmaPoints, type)
%NONLINRTSSMOOTHER Filters measurement sequence Y using a 
% non-linear Kalman filter. 
%
%Input:
%   Y           [m x N] Measurement sequence for times 1,...,N
%   x_0         [n x 1] Prior mean for time 0
%   P_0         [n x n] Prior covariance
%   f                   Motion model function handle
%   T                   Sampling time
%   Q           [n x n] Process noise covariance
%   S           [n x N] Sensor position vector sequence
%   h                   Measurement model function handle
%   R           [n x n] Measurement noise covariance
%   sigmaPoints Handle to function that generates sigma points.
%   type        String that specifies type of non-linear filter/smoother
%
%Output:
%   xf          [n x N]     Filtered estimates for times 1,...,N
%   Pf          [n x n x N] Filter error convariance
%   xp          [n x N]     Predicted estimates for times 1,...,N
%   Pp          [n x n x N] Filter error convariance
%   xs          [n x N]     Smoothed estimates for times 1,...,N
%   Ps          [n x n x N] Smoothing error convariance

% your code here!
% We have offered you functions that do the non-linear Kalman prediction and update steps.
% Call the functions using
% [xPred, PPred] = nonLinKFprediction(x_0, P_0, f, T, Q, sigmaPoints, type);
% [xf, Pf] = nonLinKFupdate(xPred, PPred, Y, S, h, R, sigmaPoints, type);

  % Size of the state and number of time steps
    n = length(x_0);
    N = size(Y, 2);

    % Allocate memory
    xs = zeros(n, N);
    Ps = zeros(n, n, N);
    xf = zeros(n, N);
    Pf = zeros(n, n, N);
    xp = zeros(n, N);
    Pp = zeros(n, n, N);
    
    % Filter
    [xp(:,1), Pp(:,:,1)] = nonLinKFprediction(x_0, P_0, f, T, Q, sigmaPoints, type);
    for k = 1:N
        if k > 1
            [xp(:,k), Pp(:,:,k)] = nonLinKFprediction(xf(:,k-1), Pf(:,:,k-1), f, T, Q, sigmaPoints, type);
        end
        [xf(:,k), Pf(:,:,k)] = nonLinKFupdate(xp(:,k), Pp(:,:,k), Y(:,k), S(:,k), h, R, sigmaPoints, type);
    end

    % Smoother
    xs(:,N) = xf(:,N);
    Ps(:,:,N) = Pf(:,:,N);
    for k = (N-1):-1:1
        [xs(:,k), Ps(:,:,k)] = nonLinRTSSupdate(xs(:,k+1), Ps(:,:,k+1), xf(:,k), Pf(:,:,k), xp(:,k+1), Pp(:,:,k+1), f, T, sigmaPoints, type);
    end


end

function [xs, Ps] = nonLinRTSSupdate(xs_kplus1, Ps_kplus1, xf_k, Pf_k, xp_kplus1, Pp_kplus1, ...
                                     f, T, sigmaPoints, type)
    % Your code here! Copy from previous task!
     n = length(xf_k);
if strcmp(type,'EKF')
    [~, F] = f(xf_k, T);
    Pxf = Pf_k * F';
else
            [SP, W] = sigmaPoints(xf_k, Pf_k, type);
            SP_kplus1 = zeros(n, size(SP, 2));
            for i = 1:size(SP, 2)
                SP_kplus1(:, i) = f(SP(:, i), T);
            end
            
            Pxf = zeros(n, n);
            for i = 1:size(SP, 2)
                Pxf = Pxf + W(i) * (SP(:, i) - xf_k) * (SP_kplus1(:, i) - xp_kplus1)';
            end
            
end

    % Compute smoother gain
    Ks = Pxf / Pp_kplus1;

    % Update the smoothed state and covariance
    xs = xf_k + Ks * (xs_kplus1 - xp_kplus1);
    Ps = Pf_k - Ks * (Pp_kplus1-Ps_kplus1) * Ks';
end

function [x, P] = nonLinKFprediction(x, P, f, T, Q, sigmaPoints, type)
%NONLINKFPREDICTION calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   f           Motion model function handle
%   T           Sampling time
%   Q           [n x n] Process noise covariance
%   sigmaPoints Handle to function that generates sigma points.
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] predicted state mean
%   P           [n x n] predicted state covariance
%

    switch type
        case 'EKF'

            % Evaluate motion model
            [fx, Fx] = f(x,T);
            % State prediction
            x = fx;
            % Covariance prediciton
            P = Fx*P*Fx' + Q;
            % Make sure P is symmetric
            P = 0.5*(P + P');

        case 'UKF'

            % Predict
            [x, P] = predictMeanAndCovWithSigmaPoints(x, P, f, T, Q, sigmaPoints, type);

            if min(eig(P))<=0
                [v,e] = eig(P);
                emin = 1e-3;
                e = diag(max(diag(e),emin));
                P = v*e*v';
            end

        case 'CKF'

            % Predict
            [x, P] = predictMeanAndCovWithSigmaPoints(x, P, f, T, Q, sigmaPoints, type);

        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end
end

function [x, P] = nonLinKFupdate(x, P, y, s, h, R, sigmaPoints, type)
%NONLINKFUPDATE calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   y           [m x 1] measurement vector
%   s           [2 x 1] sensor position vector
%   h           Measurement model function handle
%   R           [n x n] Measurement noise covariance
%   sigmaPoints Handle to function that generates sigma points.
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] updated state mean
%   P           [n x n] updated state covariance
%


switch type
    case 'EKF'
        
        % Evaluate measurement model
        [hx, Hx] = h(x,s);
        
        % Innovation covariance
        S = Hx*P*Hx' + R;
        % Kalman gain
        K = (P*Hx')/S;
        
        % State update
        x = x + K*(y - hx);
        % Covariance update
        P = P - K*S*K';
        
        % Make sure P is symmetric
        P = 0.5*(P + P');
        
    case 'UKF'

        % Update mean and covariance
        [x, P] = updateMeanAndCovWithSigmaPoints(x, P, y, s, h, R, sigmaPoints, type);
        
        if min(eig(P))<=0
            [v,e] = eig(P);
            emin = 1e-3;
            e = diag(max(diag(e),emin));
            P = v*e*v';
        end
        
    case 'CKF'

        % Update mean and covariance
        [x, P] = updateMeanAndCovWithSigmaPoints(x, P, y, s, h, R, sigmaPoints, type);
        
    otherwise
        error('Incorrect type of non-linear Kalman filter')
end

end


function [x, P] = predictMeanAndCovWithSigmaPoints(x, P, f, T, Q, sigmaPoints, type)
%
%PREDICTMEANANDCOVWITHSIGMAPOINTS computes the predicted mean and covariance
%
%Input:
%   x           [n x 1] mean vector
%   P           [n x n] covariance matrix 
%   f           measurement model function handle
%   T           sample time
%   Q           [m x m] process noise covariance matrix
%
%Output:
%   x           [n x 1] Updated mean
%   P           [n x n] Updated covariance
%

    % Compute sigma points
    [SP,W] = sigmaPoints(x, P, type);

    % Dimension of state and number of sigma points
    [n, N] = size(SP);

    % Allocate memory
    fSP = zeros(n,N);

    % Predict sigma points
    for i = 1:N
        [fSP(:,i),~] = f(SP(:,i),T);
    end

    % Compute the predicted mean
    x = sum(fSP.*repmat(W,[n, 1]),2);

    % Compute predicted covariance
    P = Q;
    for i = 1:N
        P = P + W(i)*(fSP(:,i)-x)*(fSP(:,i)-x)';
    end

    % Make sure P is symmetric
    P = 0.5*(P + P');

end

function [x, P] = updateMeanAndCovWithSigmaPoints(x, P, y, s, h, R, sigmaPoints, type)
%
%UPDATEGAUSSIANWITHSIGMAPOINTS computes the updated mean and covariance
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   y           [m x 1] measurement
%   s           [2 x 1] sensor position
%   h           measurement model function handle
%   R           [m x m] measurement noise covariance matrix
%
%Output:
%   x           [n x 1] Updated mean
%   P           [n x n] Updated covariance
%

    % Compute sigma points
    [SP,W] = sigmaPoints(x, P, type);

    % Dimension of measurement
    m = size(R,1);

    % Dimension of state and number of sigma points
    [n, N] = size(SP);

    % Predicted measurement
    yhat = zeros(m,1);
    hSP = zeros(m,N);
    for i = 1:N
        [hSP(:,i),~] = h(SP(:,i),s);
        yhat = yhat + W(i)*hSP(:,i);
    end

    % Cross covariance and innovation covariance
    Pxy = zeros(n,m);
    S = R;
    for i=1:N
        Pxy = Pxy + W(i)*(SP(:,i)-x)*(hSP(:,i)-yhat)';
        S = S + W(i)*(hSP(:,i)-yhat)*(hSP(:,i)-yhat)';
    end

    % Ensure symmetry
    S = 0.5*(S+S');

    % Updated mean
    x = x+Pxy*(S\(y-yhat));
    P = P - Pxy*(S\(Pxy'));

    % Ensure symmetry
    P = 0.5*(P+P');

end

function [Xr, Wr, j] = resampl(X, W)
%RESAMPLE Resample particles and output new particles and weights.
% resampled particles. 
%
%   if old particle vector is x, new particles x_new is computed as x(:,j)
%
% Input:
%   X   [n x N] Particles, each column is a particle.
%   W   [1 x N] Weights, corresponding to the samples
%
% Output:
%   Xr  [n x N] Resampled particles, each corresponding to some particle 
%               from old weights.
%   Wr  [1 x N] New weights for the resampled particles.
%   j   [1 x N] vector of indices refering to vector of old particles

% Your code here!
N = length(W);  % Number of particles
    W = W / sum(W);  % Normalize the weights

    % Cumulative sum of weights
    C = cumsum(W);

    % Initialize the step and the current position on the roulette wheel
    u1 = rand / N;
    i = 1;
    j = zeros(1, N);

    % Main loop
    for m = 1:N
        u_m = u1 + (m - 1) / N;  % Deterministic intervals
        while u_m > C(i)
            i = i + 1;
        end
        j(m) = i;
    end

    % Draw the samples
    Xr = X(:, j);
    Wr = repmat(1/N, 1, N);  % Weights are uniform after resampling
end


function [X_k, W_k] = pfFilterStep(X_kmin1, W_kmin1, yk, proc_f, proc_Q, meas_h, meas_R)
%PFFILTERSTEP Compute one filter step of a SIS/SIR particle filter.
%
% Input:
%   X_kmin1     [n x N] Particles for state x in time k-1
%   W_kmin1     [1 x N] Weights for state x in time k-1
%   y_k         [m x 1] Measurement vector for time k
%   proc_f      Handle for process function f(x_k-1)
%   proc_Q      [n x n] process noise covariance
%   meas_h      Handle for measurement model function h(x_k)
%   meas_R      [m x m] measurement noise covariance
%
% Output:
%   X_k         [n x N] Particles for state x in time k
%   W_k         [1 x N] Weights for state x in time k

% Your code here!
% Extract number of particles
    N = size(X_kmin1, 2);

    % Initialize new particle and weight arrays
    X_k = zeros(size(X_kmin1));
    W_k = zeros(size(W_kmin1));

    % Loop over each particle
    for i = 1:N
        % Propagation/Prediction Step
        % Sample a process noise
        v_k = mvnrnd(zeros(size(proc_Q,1),1), proc_Q)';
        
        % Apply process model to particle
        X_k(:, i) = proc_f(X_kmin1(:, i)) + v_k;

        % Update/Weighting Step
        % Compute predicted measurement
        yhat_k = meas_h(X_k(:, i));
        
        % Compute weight for this particle
        % Assuming Gaussian measurement noise
        W_k(i) = W_kmin1(i) * mvnpdf(yk, yhat_k, meas_R);
    end

    % Normalize the weights so they sum to 1
    W_k = W_k / sum(W_k);
end

function [xfp, Pfp, Xp, Wp] = pfFilter(x_0, P_0, Y, proc_f, proc_Q, meas_h, meas_R, ...
                             N, bResample, plotFunc)
%PFFILTER Filters measurements Y using the SIS or SIR algorithms and a
% state-space model.
%
% Input:
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   Y           [m x K] Measurement sequence to be filtered
%   proc_f      Handle for process function f(x_k-1)
%   proc_Q      [n x n] process noise covariance
%   meas_h      Handle for measurement model function h(x_k)
%   meas_R      [m x m] measurement noise covariance
%   N           Number of particles
%   bResample   boolean false - no resampling, true - resampling
%   plotFunc    Handle for plot function that is called when a filter
%               recursion has finished.
% Output:
%   xfp         [n x K] Posterior means of particle filter
%   Pfp         [n x n x K] Posterior error covariances of particle filter
%   Xp          [n x N x K] Non-resampled Particles for posterior state distribution in times 1:K
%   Wp          [N x K] Non-resampled weights for posterior state x in times 1:K

% Your code here, please. 
% If you want to be a bit fancy, then only store and output the particles if the function
% is called with more than 2 output arguments.
     % Initialize variables
    n = length(x_0);
    K = size(Y, 2);
    xfp = zeros(n, K);
    Pfp = zeros(n, n, K);

    % Check if we need to store particles and weights
    if nargout > 2
        Xp = zeros(n, N, K);
        Wp = zeros(N, K);
    end

    % Initialize particles and weights
    Xk = mvnrnd(x_0, P_0, N)';
    Wk = ones(1, N) / N;

    % Iterate through the measurement sequence
    for k = 1:K
        % Perform the filter step
        [Xk, Wk] = pfFilterStep(Xk, Wk, Y(:, k), proc_f, proc_Q, meas_h, meas_R);

        % Resample particles if necessary
        if bResample
            [Xk, Wk] = resampl(Xk, Wk);
        end

        % Compute the posterior mean and covariance
        xfp(:, k) = Xk * Wk';
        Pfp(:, :, k) = (Xk - xfp(:, k)) * diag(Wk) * (Xk - xfp(:, k))';

        % Store particles and weights if needed
        if nargout > 2
            Xp(:, :, k) = Xk;
            Wp(:, k) = Wk;
        end

        % Call the plot function if provided
        if ~isempty(plotFunc)
            plotFunc(k, xfp(:, k), Pfp(:, :, k), Xk, Wk);
        end
    end
end

function [fx, Fx] = coordinatedTurnMotion(x, T)


px = x(1);                  % Initialize new variables as states variables for better documentation
py = x(2);                  % Initialize new variables as states variables for better documentation
v = x(3);                   % Initialize new variables as states variables for better documentation
phi = x(4);                 % Initialize new variables as states variables for better documentation
omega = x(5);               % Initialize new variables as states variables for better documentation

% Your code for the motion model here
fx = [                      % Motion Model
      px + T*v*cos(phi);
      py + T*v*sin(phi);
      v;
      phi + T*omega;
      omega;
      ];

%Check if the Jacobian is requested by the calling function
if nargout > 1
    % Your code for the motion model Jacobian here
    Fx = [                  % Jacobian of the Motion Model
          1, 0, T*cos(phi), -T*v*sin(phi), 0;
          0, 1, T*sin(phi), T*v*cos(phi), 0;
          0, 0, 1, 0, 0;
          0, 0, 0, 1, T;
          0, 0, 0, 0, 1;
         ];
end
end
function Y = genNonLinearMeasurementSequence(X, h, R)
%GENNONLINEARMEASUREMENTSEQUENCE generates ovservations of the states 
% sequence X using a non-linear measurement model.
%
%Input:
%   X           [n x N+1] State vector sequence
%   h           Measurement model function handle
%   h           Measurement model function handle
%               [hx,Hx]=h(x) 
%               Takes as input x (state) 
%               Returns hx and Hx, measurement model and Jacobian evaluated at x
%   R           [m x m] Measurement noise covariance
%
%Output:
%   Y           [m x N] Measurement sequence
%

% Your code here
Y = zeros(size(R,1),size(X,2)-1);                       % Initialize observation sequence
r = (mvnrnd(zeros(size(R,1),1),R,size(X,2)-1))';        % Draw measurement noise from multidimensional Gaussian Distribution
for iterator = 1:(size(X,2)-1)
    [hx,~] = h(X(:,iterator+1));                        % Generate observations 
    Y(:,iterator) = hx+r(:,iterator);                   % Generate observation sequence
end
end

function [h, H] = rangeBearingMeasurements(x, s)
%RANGEBEARINGMEASUREMENTS calculates the range and the bearing to the
%position given by the state vector x, from a sensor locateed in s
%
%Input:
%   x           [n x 1] State vector
%   s           [2 x 1] Sensor position
%
%Output:
%   h           [2 x 1] measurement vector
%   H           [2 x n] measurement model Jacobian
%
% NOTE: the measurement model assumes that in the state vector x, the first
% two states are X-position and Y-position.

    % Range
    rng = norm(x(1:2)-s);
    % Bearing
    ber = atan2(x(2)-s(2),x(1)-s(1));
    % Measurement vector
    h = [rng;ber];

    % Measurement model Jacobian
    H = [
        (x(1)-s(1))/rng      (x(2)-s(2))/rng     0 0 0;
        -(x(2)-s(2))/(rng^2) (x(1)-s(1))/(rng^2) 0 0 0
        ];

end

function [x, y] = samplesToState(measurements, S)
    distance = measurements(1,:);
    angle = measurements(2,:);
    
    pos = S + distance .* [cos(angle); sin(angle)];
    x = pos(1,:);
    y = pos(2,:);
    
end

function [SP,W] = sigmaPoints(x, P, type)
% SIGMAPOINTS computes sigma points, either using unscented transform or
% using cubature.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%
%Output:
%   SP          [n x 2n+1] UKF, [n x 2n] CKF. Matrix with sigma points
%   W           [1 x 2n+1] UKF, [1 x 2n] UKF. Vector with sigma point weights 
%
n = size(x,1);
P_sqrt = sqrtm(P);

switch type        
    case 'UKF'
        num_points = 2*n + 1;
        
        SP = zeros(n, num_points);
        SP(:,1) = x;
        
        W = zeros(1, num_points);
        W0 = 1 - n/3;
        
        W(1) = W0;
        SP(:,1) = x;
        
        for i = 1:n
           W(i+1) = (1-W0) / (2*n); 
           W(i+1+n) = (1-W0) / (2*n); 
           
           SP(:,i+1)   = x + sqrt(n/(1-W0))*P_sqrt(:,i);
           SP(:,i+1+n) = x - sqrt(n/(1-W0))*P_sqrt(:,i);
        end
        
            
    case 'CKF'
        num_points = 2*n;
        
        SP = zeros(n, num_points);
        W = zeros(1, num_points);
        
        for i = 1:n
            W(i)      = 1/(2*n);
            W(i+n)    = W(i); 
            SP(:,i)   = x + sqrt(n)*P_sqrt(:,i);
            SP(:,i+n) = x - sqrt(n)*P_sqrt(:,i);
        end
        
    otherwise
        error('Incorrect type of sigma point')
end

end


function X = genLinearStateSequence(x_0, P_0, A, Q, N)
%GENLINEARSTATESEQUENCE generates an N-long sequence of states using a 
%    Gaussian prior and a linear Gaussian process model
%
%Input:
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   A           [n x n] State transition matrix
%   Q           [n x n] Process noise covariance
%   N           [1 x 1] Number of states to generate
%
%Output:
%   X           [n x N+1] State vector sequence
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Your code here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate initial state from prior.
X0 = mvnrnd(x_0, P_0, 1)';

% Generate state sequence by letting the initial state propagate
X = zeros(length(x_0), N+1);
X(:,1) = X0;

for k = 2:N+1
    X(:,k) = mvnrnd(A * X(:,k-1), Q, 1);
end


end

function [X, P, V] = kalmanFilter(Y, x_0, P_0, A, Q, H, R)
%KALMANFILTER Filters measurements sequence Y using a Kalman filter. 
%
%Input:
%   Y           [m x N] Measurement sequence
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   A           [n x n] State transition matrix
%   Q           [n x n] Process noise covariance
%   H           [m x n] Measurement model matrix
%   R           [m x m] Measurement noise covariance
%
%Output:
%   x           [n x N] Estimated state vector sequence
%   P           [n x n x N] Filter error convariance
%   V           [m x N] inovation

% Parameters
N = size(Y,2);

n = length(x_0);
m = size(Y,1);

% Data allocation
X = zeros(n,N);
P = zeros(n,n,N);
V = zeros(1,N);

% 1. Predict the next state
% 2. Update the prediction with measurement
% 3. Save the update for next iteration

for k = 1:N
    [x_0, P_0] = linearPrediction(x_0, P_0, A, Q);
    [x_0, P_0, V_0] = linearUpdate(x_0, P_0, Y(:,k), H, R);
    V(:, k) = V_0;
    X(:,k) = x_0;
    P(:,:,k) = P_0;
end
end

function [x, P] = linearPrediction(x, P, A, Q)
%LINEARPREDICTION calculates mean and covariance of predicted state
%   density using a liear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   A           [n x n] State transition matrix
%   Q           [n x n] Process noise covariance
%
%Output:
%   x           [n x 1] predicted state mean
%   P           [n x n] predicted state covariance
%

% Your code here

% Calculate posterior from the process model
x = A*x;
P = A*P*A' + Q;

end


function [x, P, V] = linearUpdate(x, P, y, H, R)
%LINEARPREDICTION calculates mean and covariance of predicted state
%   density using a linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   y           [m x 1] Measurement
%   H           [m x n] Measurement model matrix
%   R           [m x m] Measurement noise covariance
%
%Output:
%   x           [n x 1] updated state mean
%   P           [n x n] updated state covariance
%   V           [m x 1] inovation

% Your code here
S = H*P*H' + R;
K = P*H'*inv(S);
V = y - H*x;

x = x + K*V;
P = P - K*S*K';

end






































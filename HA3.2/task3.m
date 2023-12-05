%% task 3: Tuning non-linear filters
clc;
clear;
% a b: Tune Q and evaluate the performance of the filter
sigma_v = 0.001 * 1e-4;
sigma_w = 0.001*pi/180;

% True track
% Sampling period
T = 0.1;
% Length of time sequence
K = 600;
% Allocate memory
omega = zeros(1,K+1);
% Turn rate
omega(200:400) = -pi/201/T;
% Initial state
x0 = [0 0 20 0 omega(1)]';
% Allocate memory
X = zeros(length(x0),K+1);
X(:,1) = x0;
% Create true track
for i=2:K+1
    % Simulate
    X(:,i) = coordinatedTurnMotion(X(:,i-1), T);
    % Set turn-rate
    X(5,i) = omega(i);
end

x_0 = [0,0,0,0,0]';
P_0 = diag([10,10,10,(5*pi/180),(1*pi/180)].^2);

s1 = [300,-100]';
s2 = [300,-300]';

sigma_phi1 = pi/180;
sigma_phi2 = pi/180;

R = diag([sigma_phi1, sigma_phi2].^2);
h = @(x) dualBearingMeasurement(x,s1,s2);
Y = genNonLinearMeasurementSequence(X,h,R);

f = @(x) coordinatedTurnMotion(x,T);
Q = diag([0 0 T*sigma_v^2 0 T*sigma_w^2]);

[xf, Pf, xp, Pp] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q, h, R, 'CKF');

% Calcualte unfiltered position from sensors given angles
Xm(1,:) = (s2(2)-s1(2)+tan(Y(1,:))*s1(1)-tan(Y(2,:))*s2(1))./(tan(Y(1,:))- tan(Y(2,:)));
Xm(2,:) = s1(2)+tan(Y(1,:)).*(Xm(1,:)-s1(1));


% c
% First figure
figure();
clf;
hold on
grid on
plot(s1(1), s1(2), '*c', 'Linewidth', 2);
plot(s2(1), s2(2), 'oy', 'Linewidth', 2);
plot(X(1,:), X(2,:), '-r', 'Linewidth', 1.5);
plot(xf(1,:), xf(2,:), ':g', 'Linewidth', 1.5);
plot(Xm(1,:), Xm(2,:), 'xb', 'MarkerSize', 3);
for i = 1:10:length(xf)
    [xy4] = sigmaEllipse2D(xf(1:2, i), Pf(1:2, 1:2, i), 3, 50);
    plot(xy4(1,:), xy4(2,:), '--m', 'Linewidth', 1);
end
xlabel('x');
ylabel('y');
legend('Sensor 1 Position', 'Sensor 2 Position', 'True State', 'Filtered State', 'Measured State', '$3\sigma$ region', 'Interpreter', 'Latex', 'Location', 'Southeast');
title("Tuning of Kalman Filter, Low Process Noise Case");
axis equal
hold off

% Second figure
figure();
clf;
hold on
grid on
plot((1:K) * T, vecnorm(xf(1:2, :) - X(1:2, 2:end), 2, 1), 'LineWidth', 1.5, 'Color', [0.5, 0, 0.5]);
xlabel('Time [s]');
ylabel('$||p_{k}-\hat{p}_{k|k}||_{2}$', 'Interpreter', 'Latex');
title("Tuning of Kalman Filter, Low Process Noise Case");
hold off

















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
    switch type        
        case 'UKF'
    
            SP = zeros(n,2*n+1);                                    % Initialize Sigma Points            
            W0 = 1-n/3;                                             % Initialize Weights
            SP(:,1) = x;                                            % Initialize First Sigma Point
            SqrtP = sqrtm(P);
            
            for i = 1:size(P,2)
                SP(:,i+1) = x + sqrt((n)/(1-W0))*SqrtP(:,i);        % Calculate Remaining Sigma Points
                SP(:,i+n+1) = x - sqrt((n)/(1-W0))*SqrtP(:,i);      % Calculate Remaining Sigma Points
            end
            
            W = [W0, ((1-W0)/(2*n))*ones(1,2*n)];                   % Calculate Remaining Weights
                
        case 'CKF'
            
            SP = zeros(n,2*n);                                      % Initialize Sigma Points
            SqrtP = sqrtm(P);
            
            for i = 1:size(P,2)
                SP(:,i) = x + sqrt(n)*SqrtP(:,i);                   % Calculate Sigma Points
                SP(:,i+n) = x - sqrt(n)*SqrtP(:,i);                 % Calculate Sigma Points
            end
            
            W = ((1)/(2*n))*ones(1,2*n);                            % Calculate Weights
            
        otherwise
            error('Incorrect type of sigma point')
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

function [hx, Hx] = dualBearingMeasurement(x, s1, s2)
%DUOBEARINGMEASUREMENT calculates the bearings from two sensors, located in 
%s1 and s2, to the position given by the state vector x. Also returns the
%Jacobian of the model at x.
%
%Input:
%   x           [n x 1] State vector, the two first element are 2D position
%   s1          [2 x 1] Sensor position (2D) for sensor 1
%   s2          [2 x 1] Sensor position (2D) for sensor 2
%
%Output:
%   hx          [2 x 1] measurement vector
%   Hx          [2 x n] measurement model Jacobian
%
% NOTE: the measurement model assumes that in the state vector x, the first
% two states are X-position and Y-position.

% Your code here
n = size(x,1);      
N = size(x,2);
hx = zeros(2,N);        % Initialize measurement vector with zeros
Hx = zeros(2,n);        % Initialize measurement model Jacobian with zeros

hx(1:2,:) = [                  % Measurement vector
     atan2(x(2,:)-s1(2),x(1,:)-s1(1));
     atan2(x(2,:)-s2(2),x(1,:)-s2(1));
	 ];
        
Hx(1:2,1:2) = [         % Measurement Model Jacobian
     -(x(2)-s1(2))/( (x(1)-s1(1))^2 + (x(2)-s1(2))^2 ), (x(1)-s1(1))/( (x(1)-s1(1))^2 + (x(2)-s1(2))^2 );
     -(x(2)-s2(2))/( (x(1)-s2(1))^2 + (x(2)-s2(2))^2 ), (x(1)-s2(1))/( (x(1)-s2(1))^2 + (x(2)-s2(2))^2 );
     ];
 
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
function [xf, Pf, xp, Pp] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q, h, R, type)
    % NONLINEARKALMANFILTER Filters measurement sequence Y using a 
    % non-linear Kalman filter. 
    %
    % Input:
    %   Y           [m x N] Measurement sequence for times 1,...,N
    %   x_0         [n x 1] Prior mean for time 0
    %   P_0         [n x n] Prior covariance
    %   f                   Motion model function handle
    %                       [fx,Fx]=f(x) 
    %                       Takes as input x (state) 
    %                       Returns fx and Fx, motion model and Jacobian evaluated at x
    %   Q           [n x n] Process noise covariance
    %   h                   Measurement model function handle
    %                       [hx,Hx]=h(x,T) 
    %                       Takes as input x (state), 
    %                       Returns hx and Hx, measurement model and Jacobian evaluated at x
    %   R           [m x m] Measurement noise covariance
    %
    % Output:
    %   xf          [n x N]     Filtered estimates for times 1,...,N
    %   Pf          [n x n x N] Filter error convariance
    %   xp          [n x N]     Predicted estimates for times 1,...,N
    %   Pp          [n x n x N] Predicted error convariance

% Parameters
N = size(Y,2);
n = length(x_0);
m = size(Y,1);

% Data allocation
xp = zeros(n,N);
Pp = zeros(n,n,N);

xf = zeros(n,N+1);
Pf = zeros(n,n,N+1);

% Filter Implementation

xf(:,1)   = x_0;
Pf(:,:,1) = P_0;

for iterator = 1:N
    [xp(:,iterator), Pp(:,:,iterator)] = nonLinKFprediction(xf(:,iterator), Pf(:,:,iterator), f, Q, type);
    [xf(:,iterator+1), Pf(:,:,iterator+1)] = nonLinKFupdate(xp(:,iterator), Pp(:,:,iterator), Y(:,iterator), h, R, type);
end
    
xf = xf(:,2:end);
Pf = Pf(:,:,2:end);
    
end

function [x, P] = nonLinKFprediction(x, P, f, Q, type)
%NONLINKFPREDICTION calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   f           Motion model function handle
%               [fx,Fx]=f(x) 
%               Takes as input x (state), 
%               Returns fx and Fx, motion model and Jacobian evaluated at x
%               All other model parameters, such as sample time T,
%               must be included in the function
%   Q           [n x n] Process noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] predicted state mean
%   P           [n x n] predicted state covariance
%

    n = size(x,1);
    switch type
        case 'EKF'
            
            % Your EKF code here
            [fx,Fx] = f(x);                 
            x = fx;                                                 % Calculate Predicted Mean
            P = Fx*P*Fx' + Q;                                       % Calculate Predicted Covariance
            
        case 'UKF'
    
            % Your UKF code here
            [SP,W] = sigmaPoints(x, P, type);                       % Calculate Sigma Points
            x = zeros(n,1);                                         % Initialize Predicted Mean
            for i = 1:size(W,2)
                x = x + f(SP(:,i))*W(i);                            % Calculate Predicted Mean
            end
            
            P = Q;                                                  % Initialize Predicted Covariance
            for i = 1:size(W,2)
                P = P + (f(SP(:,i)) - x)*(f(SP(:,i)) - x).'*W(i);   % Calculate Predicted Covariance
            end
            
            % Make sure the covariance matrix is semi-definite
            if min(eig(P))<=0
                [v,e] = eig(P, 'vector');
                e(e<0) = 1e-4;
                P = v*diag(e)/v;
            end
                
        case 'CKF'
            
            % Your CKF code here
            [SP,W] = sigmaPoints(x, P, type);                       % Calculate Sigma Points
            x = zeros(n,1);                                         % Initialize Predicted Mean
            for i = 1:size(W,2)
                x = x + f(SP(:,i))*W(i);                            % Calculate Predicted Mean
            end
            
            P = Q;                                                  % Initialize Predicted Covariance
            for i = 1:size(W,2)
                P = P + (f(SP(:,i)) - x)*(f(SP(:,i)) - x).'*W(i);   % Calculate Predicted Covariance
            end
            
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end

end



function [x, P] = nonLinKFupdate(x, P, y, h, R, type)
%NONLINKFUPDATE calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   y           [m x 1] measurement vector
%   h           Measurement model function handle
%               [hx,Hx]=h(x) 
%               Takes as input x (state), 
%               Returns hx and Hx, measurement model and Jacobian evaluated at x
%               Function must include all model parameters for the particular model, 
%               such as sensor position for some models.
%   R           [m x m] Measurement noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] updated state mean
%   P           [n x n] updated state covariance
%

    switch type
        case 'EKF'
            
            % Your EKF update here
            [hx,Hx] = h(x);            
            S = Hx*P*Hx' + R;                                           % Calculate Innovation Covariance
            K = P*Hx'*inv(S);                                           % Calculate Innovation Gain
            x = x + K*(y - hx);                                         % Updated Mean
            P = P - K*S*K';                                             % Updated Covariance
        
        case 'UKF'
    
            % Your UKF update here
            [SP,W] = sigmaPoints(x, P, type);
            yhat = zeros(size(y));            
            for i = 1:size(W,2)
                yhat = yhat + h(SP(:,i))*W(i);                          % Calculate Measurements
            end
            
            Pxy = zeros(size(P,1),size(R,1));            
            for i = 1:size(W,2)
                Pxy = Pxy + (SP(:,i) - x)*(h(SP(:,i)) - yhat)'*W(i);    % Calculate Cross-Covariance
            end

            S = R;
            for i = 1:size(W,2)
                S = S + (h(SP(:,i)) - yhat)*(h(SP(:,i)) - yhat)'*W(i);  % Calculate Innovation Covariance
            end
            
            x = x + Pxy*inv(S)*(y - yhat);                              % Updated Mean
            P = P - Pxy*inv(S)*Pxy';                                    % Updated Covariance
    
            % Make sure the covariance matrix is semi-definite
            if min(eig(P))<=0
                [v,e] = eig(P, 'vector');
                e(e<0) = 1e-4;
                P = v*diag(e)/v;
            end
            
        case 'CKF'
    
            % Your CKF update here
            [SP,W] = sigmaPoints(x, P, type);
            yhat = zeros(size(y));            
            for i = 1:size(W,2)
                yhat = yhat + h(SP(:,i))*W(i);                          % Calculate Measurements
            end
            
            Pxy = zeros(size(P,1),size(R,1));            
            for i = 1:size(W,2)
                Pxy = Pxy + (SP(:,i) - x)*(h(SP(:,i)) - yhat)'*W(i);    % Calculate Cross-Covariance
            end

            S = R;
            for i = 1:size(W,2)
                S = S + (h(SP(:,i)) - yhat)*(h(SP(:,i)) - yhat)'*W(i);  % Calculate Innovation Covariance
            end
            
            x = x + Pxy*inv(S)*(y - yhat);                              % Updated Mean
            P = P - Pxy*inv(S)*Pxy';                                    % Updated Covariance
    
            % Make sure the covariance matrix is semi-definite
            if min(eig(P))<=0
                [v,e] = eig(P, 'vector');
                e(e<0) = 1e-4;
                P = v*diag(e)/v;
            end
            
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end

end


function [ xy ] = sigmaEllipse2D( mu, Sigma, level, npoints )
%SIGMAELLIPSE2D generates x,y-points which lie on the ellipse describing
% a sigma level in the Gaussian density defined by mean and covariance.
%
%Input:
%   MU          [2 x 1] Mean of the Gaussian density
%   SIGMA       [2 x 2] Covariance matrix of the Gaussian density
%   LEVEL       Which sigma level curve to plot. Can take any positive value, 
%               but common choices are 1, 2 or 3. Default = 3.
%   NPOINTS     Number of points on the ellipse to generate. Default = 32.
%
%Output:
%   XY          [2 x npoints] matrix. First row holds x-coordinates, second
%               row holds the y-coordinates. First and last columns should 
%               be the same point, to create a closed curve.


%Setting default values, in case only mu and Sigma are specified.
if nargin < 3
    level = 3;
end
if nargin < 4
    npoints = 32;
end

phi = linspace(0,2*pi,npoints);	% Calculate an even spaced vector of angles
z = level*[cos(phi);sin(phi)];	% Calculate the new vector points for the ellipse
xy = mu + sqrtm(Sigma)*z;	% Generate 2D points lying on a specific level curve given density
end


















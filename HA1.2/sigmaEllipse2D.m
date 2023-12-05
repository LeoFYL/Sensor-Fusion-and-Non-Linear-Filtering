% function sigmaEllipse2D(mu, Sigma, nSigma)
% %SIGHMAELLIPSE2D plots the 2D sigma ellipse of a Gaussian distribution
% %   SIGHMAELLIPSE2D(MU, SIGMA, NSIGMA) plots the 2D sigma ellipse of a
% %   Gaussian distribution specified by its mean MU and covariance SIGMA.
% %   The ellipse is drawn with a radius of NSIGMA standard deviations.
% 
% 
% % Generate points on unit circle
% theta = linspace(0, 2*pi, 100)';
% xy = [cos(theta), sin(theta)];
% 
% % Calculate eigenvectors and eigenvalues of Sigma
% [V, D] = eig(Sigma);
% 
% % Calculate points on ellipse
% xy = nSigma * (V * sqrt(D)) * xy';
% 
% % Plot ellipse
% plot(mu(1) + xy(1,:), mu(2) + xy(2,:), 'LineWidth', 2);
% 
% 
% 
% end

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

    %Your code here

    %Evenly spaced points
    theta = linspace(0, 2*pi, npoints);
    
    
    % Make sure the covariance matrix is semi-definite
    if min(eig(Sigma))<=0
        [v,e] = eig(Sigma, 'vector');
        e(e<0) = 1e-4;
        Sigma = v*diag(e)/v;
    end

    %Level curve
    xy = mu + level*sqrtm(Sigma) * [cos(theta); sin(theta)];

end


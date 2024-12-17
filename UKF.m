function [xhat, P, innovation, S] = UKF(xhat_prev, P_prev, u, y, dt, Q, R, alpha, beta, kappa, n)
%UKF Perform a single UKF update step
% Inputs:
%   xhat_prev -> Previous state estimate
%   P_prev -> Previous state covariance
%   u -> Control input
%   y -> Measurement
%   dt -> Time step
%   Q -> Process noise covariance
%   R -> Measurement noise covariance
%   gamma -> Scaling factor for sigma points
%   lambda -> Composite scaling parameter
%   n -> State dimension
% Outputs:
%   xhat -> Updated state estimate
%   P -> Updated state covariance


% Generate Sigma Points
lambda = alpha^2 * (n + kappa) - n;
lambda2 = sqrt(n + lambda);
sigma_points = sigmaPoints(xhat_prev, P_prev, lambda2);

% Time Update (Prediction Step)
x_sigma_pred = zeros(n, 2*n + 1);
for i = 1:2*n + 1
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    [~, xpred] = ode45(@(t, x)coopEOM(t, x, u, zeros(n, 1)), [0 dt], sigma_points(:, i), options);
    x_sigma_pred(:, i) = xpred(end, :)';
end

% Predicted state mean
weights_mean = [lambda / (n + lambda), repmat(1 / (2 * (n + lambda)), 1, 2 * n)];
xhat_pred = x_sigma_pred * weights_mean';

% Predicted state covariance
weights_cov = weights_mean;
weights_cov(1) = weights_cov(1) + (1 - alpha^2 + beta);
P_pred = Q;
for i = 1:2 * n + 1
    diff1 = x_sigma_pred(:, i) - xhat_pred;
    diff1(3) = mod(diff1(3) + pi, 2*pi) - pi;
    diff1(6) = mod(diff1(6) + pi, 2*pi) - pi;

    P_pred = P_pred + weights_cov(i) * (diff1 * diff1');
end

% Measurement Update
y_sigma_pred = zeros(length(y), 2 * n + 1);
for i = 1:2 * n + 1
    y_sigma_pred(:, i) = sensors(x_sigma_pred(:, i));
end

% Predicted measurement mean
yhat_pred = y_sigma_pred * weights_mean';

% Measurement covariance
P_yy = R;
for i = 1:2 * n + 1
    diff2 = y_sigma_pred(:, i) - yhat_pred;
    diff2(1) = mod(diff2(1) + pi, 2*pi) - pi;
    diff2(3) = mod(diff2(3) + pi, 2*pi) - pi;

    P_yy = P_yy + weights_cov(i) * (diff2 * diff2');
end
S = P_yy;

% Cross-covariance
P_xy = zeros(n, length(y));
for i = 1:2 * n + 1
    diff1 = x_sigma_pred(:, i) - xhat_pred;
    diff1(3) = mod(diff1(3) + pi, 2*pi) - pi;
    diff1(6) = mod(diff1(6) + pi, 2*pi) - pi;

    diff2 = y_sigma_pred(:, i) - yhat_pred;
    diff2(1) = mod(diff2(1) + pi, 2*pi) - pi;
    diff2(3) = mod(diff2(3) + pi, 2*pi) - pi;

    P_xy = P_xy + weights_cov(i) * diff1 * diff2';
end

% Kalman gain
K = P_xy / P_yy;

% Update state and covariance
innovation = y - yhat_pred;
innovation(1) = mod(innovation(1) + pi, 2*pi) - pi;
innovation(3) = mod(innovation(3) + pi, 2*pi) - pi;

xhat = xhat_pred + K * innovation;
P = P_pred - K * P_yy * K';
end
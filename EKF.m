function [xhat_meas, P_meas] = EKF(xhat_prev, P_prev, u, y, dt, Q, R)
%EKF This performs a single EKF update with prediction and measurement
%update steps
%
% Inputs: 
%   xhat_prev -> previous estimate of the full state
%   P_prev -> previous covariance
%   u -> control input
%   y -> measurement vector
%   dt -> time step
%   Q -> process noise covariance matrix
%   R -> measurement noise covariance matrix
%
% Outputs:
%   xhat_meas -> next full state estimate
%   P_meas -> next covariance
%
% Author: Thomas Dunnington
% Modified: 12/3/2024

%% Time update step
% Deterministic state prediction
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[~, xpred] = ode45(@(t, x)coopEOM(t, x, u, zeros(length(xhat_prev),1)), [0 dt], xhat_prev, options);
xhat_pred = xpred(end, :)';

% Linearization matrices
[A, B, C, D] = linearize(xhat_prev, u); % Linearize along the current best state estimate, not the nominal trajectory
[F, G, Omega, H] = eulerDiscretize(A, B, C, D, eye(size(A)), dt);

% Covariance prediction
P_pred = F*P_prev*F' + Omega*Q*Omega';  

%% Measurement update step
% Deterministic nonlinear measurement model evaluation
yhat_pred = sensors(xhat_pred);

% Linearization matrices at new best estimate
[A, B, C, D] = linearize(xhat_pred, u);
[F, G, Omega, H] = eulerDiscretize(A, B, C, D, eye(size(A)), dt);

% Innovation vector
innovation = y - yhat_pred;

% Kalman gain
K = P_pred * H' * inv(H * P_pred * H' + R);

% Update state and covariance
xhat_meas = xhat_pred + K*innovation;
P_meas = (eye(size(F)) - K*H) * P_pred;

end


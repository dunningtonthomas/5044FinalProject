%% Main Function for EKF Analysis
close all; clear; clc;

%% EKF
Data = load('cooplocalization_finalproj_KFdata.mat');

% Ode values
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

% Test the EKF update
% Nominal values
x_ugv = [10; 0; pi/2];
x_uav = [-60; 0; -pi/2];
u_ugv = [2; -pi/18];
u_uav = [12; pi/25];

x_nom = [x_ugv; x_uav];
u_nom = [u_ugv; u_uav];

% Time step
dt = 0.1;

xhat_prev = x_nom;
u = u_nom;
y = sensors(xhat_prev);
Q = Data.Qtrue;
R = Data.Rtrue;
P_prev = [10, 0, 0, 0, 0, 0;
            0, 10, 0, 0, 0, 0;
            0, 0, pi, 0, 0, 0;
            0, 0, 0, 10, 0, 0;
            0, 0, 0, 0, 10, 0;
            0, 0, 0, 0, 0, pi];

%[xhat_meas, P_meas] = EKF(xhat_prev, P_prev, u, y, dt, Q, R);

% Truth model testing, get true state and data
n_ind = 1000;
[time_tmt, x_noise_mat, y_noise_mat] = simulateNoise(x_nom, u_nom, Q, R, dt, n_ind);

% Simulate nonlinear equations to get deterministic state trajectory


% Nonlinear least squares warm start
num_measurements = 10;

% Create sensor noise matrix
R_cells = repmat({R}, 1, num_measurements);
Rbig = blkdiag(R_cells{:});

% Get measurements of data and stack the data
y_meas = [];
for i = 1:num_measurements
    % y_meas = [y_meas; Data.ydata(:,i+1)];
    y_meas = [y_meas; y_noise_mat(i,:)'];
end

% Initial guess
xhat_old = [0; 0; 0; y_meas(4); y_meas(5); 0]; % Rough initial guess, maybe use the first measurement
%xhat_old = x_nom;

% Gauss Newton optimization
max_iterations = 200;
alpha_splits = 10;
iteration = 0;
dt = 0.1;
time_vector = 0:dt:(num_measurements*dt);

while iteration < max_iterations
    % Get predicted measurement vector
    y_pred = [];
    Hbig = [];
    xinit = xhat_old;
    [~, xpred] = ode45(@(t, x)coopEOM(t, x, u, zeros(length(xinit),1)), time_vector, xinit, options);

    for i = 1:num_measurements
        % Get predicted measurement
        y_pred = [y_pred; sensors(xpred(i+1, :)')];

        % Measurement jacobian
        [~, ~, H, ~] = linearize(xpred(i+1, :)', zeros(4,1));
        Hbig = [Hbig; H];
    end

    % Calculate the cost function
    Jcurr = (y_meas - y_pred)' * inv(Rbig) * (y_meas - y_pred);
    dxhat = inv(Hbig' * inv(Rbig) * Hbig) * Hbig' * inv(Rbig) * (y_meas - y_pred);

    % No alpha splitting
    % xhat_new = xhat_old + dxhat;

    % Make alpha smaller so we don't diverge
    alpha = 1;
    alpha_iter = 1;
    while alpha_iter < alpha_splits
        % Proposed update to the estimate
        xhat_new = xhat_old + alpha .* dxhat;

        % Get new predicted measurements
        xinit = xhat_new;
        [~, xpred] = ode45(@(t, x)coopEOM(t, x, u, zeros(length(xinit),1)), time_vector, xinit, options);

        yhat_new = [];
        for i = 1:num_measurements
            % Predicted sensor measurements
            yhat_new = [yhat_new; sensors(xpred(i+1, :)')];
        end

        % Calculate the cost function
        Jnew = (y_meas - yhat_new)' * inv(Rbig) * (y_meas - yhat_new);

        % Take best cost
        if Jnew > Jcurr
            alpha = alpha / 2; % Shrink size
        else
            break;
        end

    end

    % Output the iteration and the cost function
    fprintf('Iteration: %d, Cost: %.6f\n', iteration, Jcurr);

    % Update the best guess
    xhat_old = xhat_new;
    iteration = iteration + 1;
end

% Final estimate
xhat_final = xhat_old;

% Calculate the final covariance
[~, ~, H, ~] = linearize(xhat_final, zeros(4,1));
P_final = inv(H'*H);


% Output the final xhat measurement
disp('Final Initial State Estimate:');
disp(xhat_new);

disp('Final Covariance Estimate:');
disp(P_final);



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
num_measurements = 100;
%[xhat_final, P_final] = NLLS(u, y_noise_mat, num_measurements, Q, R, dt);

% Run the EKF with an initial guess and covariance
x_init = x_noise_mat(1,:);
%x_init = xhat_final;
P_init = diag([1000, 1000, 2*pi, 1000, 1000, 2*pi]);

% Tune the Q and R matrices
Q = 1000 .* Q;
R = R;

% Estimated and covariance matrices
xhat_mat = zeros(length(x_noise_mat(:,1)), 6);
sigma_mat = zeros(length(y_noise_mat(:,1)), 6);

% Monte Carlo Simulation Parameters
Nsim = 50; % Number of Monte Carlo simulations
Nstate = size(x_nom, 1);
Nmeas = size(y_noise_mat, 2);

% Preallocate for NEES and NIS results
nees_vals = zeros(Nsim, length(time_tmt)-1);
nis_vals = zeros(Nsim, length(time_tmt)-1);

for sim_idx = 1:Nsim
    % Initialize the EKF for each Monte Carlo run
    xhat_curr = x_init;
    P_curr = P_init;

    % Loop over the noisy measurements to get the estimates
    for i = 1:length(y_noise_mat(:,1))
        % Get the current measurement
        y_meas = y_noise_mat(i,:)';

        % Linearize the measurement function to get the Jacobian H
        [~, ~, H, ~] = linearize(xhat_curr, u_nom);

        % Kalman update
        [xhat_curr, P_curr] = EKF(xhat_curr, P_curr, u_nom, y_meas, dt, Q, R);

        % Store the result in a matrix of the estimated state
        xhat_mat(i+1, :) = xhat_curr';
        sigma_mat(i+1,:) = sqrt(diag(P_curr))';

        % Compute the NEES and NIS for this simulation
        % NEES (Normalized Estimation Error Squared)
        x_error = xhat_curr - x_noise_mat(i,:)';
        nees_vals(sim_idx, i) = x_error' * inv(P_curr) * x_error;

        % Innovation vector and covariance
        innovation = y_meas - H * xhat_curr;
        S = R + H * P_curr * H'; % Innovation covariance
        nis_vals(sim_idx, i) = innovation' * inv(S) * innovation;
    end
end

% Chi-Square Test Bounds
alpha = 0.05; % Significance level
r1_NEES = chi2inv(alpha/2, Nstate) / Nsim;
r2_NEES = chi2inv(1-alpha/2, Nstate) / Nsim;

r1_NIS = chi2inv(alpha/2, Nmeas) / Nsim;
r2_NIS = chi2inv(1-alpha/2, Nmeas) / Nsim;

% Average NEES and NIS across all runs
mean_nees = mean(nees_vals, 1); % Time-averaged NEES
mean_nis = mean(nis_vals, 1);   % Time-averaged NIS


%% Plot the results
% figure();
plotSim(time_tmt, x_noise_mat, y_noise_mat, '-');
%plotSim(time_tmt, xhat_mat, y_noise_mat, '--');


% Plot the errors with the two sigma bounds
x_error = xhat_mat - x_noise_mat;
% x_error(:,3) = mod(x_error(:,3) + pi, 2*pi) - pi;
% x_error(:,6) = mod(x_error(:,6) + pi, 2*pi) - pi;

% State labels
state_labels = {'$\xi_g$ Error [m]', '$\eta_g$ Error [m]', '$\theta_g$ Error [rad]', ...
                '$\xi_a$ Error [m]', '$\eta_a$ Error [m]', '$\theta_a$ Error [rad]'};

% Create figure
figure('Name', 'State Estimator Error', 'NumberTitle', 'off', 'Color', 'w');
plot_num = 1;

% Loop to create 3x2 subplot
plot_locations = [1, 3, 5, 2, 4, 6];
for i = 1:6
    subplot(3, 2, plot_locations(plot_num)); % 3 rows, 2 columns
    hold on;
    
    % Plot error and bounds
    plot(time_tmt(2:end), x_error(2:end,i), 'b', 'LineWidth', 1.5, 'DisplayName', 'Error'); 
    plot(time_tmt(4:end), 2*sigma_mat(4:end, i), 'r--', 'LineWidth', 1.2, 'DisplayName', '+2\sigma'); 
    plot(time_tmt(4:end), -2*sigma_mat(4:end, i), 'r--', 'LineWidth', 1.2, 'DisplayName', '-2\sigma');
    
    % Add labels and grid
    xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 10);
    ylabel(state_labels{plot_num}, 'Interpreter', 'latex', 'FontSize', 10);

    % Add legend for the top right
    if(plot_num == 4)
        legend('Estimate', '$2 \sigma$ bounds', 'interpreter', 'latex');
    end
    
    % Add title
    title([state_labels{plot_num}, ' with $\pm 2 \sigma$ Bounds'], 'Interpreter', 'latex', 'FontSize', 10);
    
    % Increment plot number
    plot_num = plot_num + 1;
end

% Add a global title
sgtitle('State Estimator Errors with $\pm 2 \sigma$  Bounds', 'Interpreter', 'latex', 'FontSize', 14);


figure;
subplot(2, 1, 1);
plot(time_tmt(2:end), mean_nees, 'b', 'LineWidth', 1.5);
hold on;
yline(r1_NEES, 'r--', 'LineWidth', 1.2);
yline(r2_NEES, 'r--', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('NEES');
legend('Mean NEES', '\chi^2 Lower Bound', '\chi^2 Upper Bound');
title('NEES Chi-Square Test');
grid on;

subplot(2, 1, 2);
plot(time_tmt(2:end), mean_nis, 'b', 'LineWidth', 1.5);
hold on;
yline(r1_NIS, 'r--', 'LineWidth', 1.2);
yline(r2_NIS, 'r--', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('NIS');
legend('Mean NIS', '\chi^2 Lower Bound', '\chi^2 Upper Bound');
title('NIS Chi-Square Test');
grid on;

sgtitle('NEES and NIS Validation for Extended Kalman Filter');




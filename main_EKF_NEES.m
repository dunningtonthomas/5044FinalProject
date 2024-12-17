%% EKF
close all; clear; clc
Data = load('cooplocalization_finalproj_KFdata.mat');

% Ode values
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

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

% TMT
n_ind = 1000;
[time_tmt, x_noise_mat, y_noise_mat] = simulateNoise(x_nom, u_nom, Q, R, dt, n_ind);

% Run the EKF with an initial guess and covariance
x_init = x_nom';
%x_init = xhat_final;
%P_init = diag([100, 100, 2*pi, 100, 100, 2*pi]);
%P_init = diag([100, 100, pi, 100, 100, pi]);
P_init = diag([10, 10, pi/4, 25, 25, pi/4]);
% P_init = diag([1000, 1000, 2*pi, 1000, 1000, 2*pi]);

% Tune the Q matrix
% Q_EKF = 10.0 .* Q;
% Q_EKF(4,4) = 1000 .* Q_EKF(4,4);
% Q_EKF(5,5) = 1000 .* Q_EKF(5,5);
Q_EKF = 1 .*Q;
Q_EKF(1,1) = 10 .* Q_EKF(1,1);
Q_EKF(2,2) = 10 .* Q_EKF(2,2);
Q_EKF(3,3) = 1 .* Q_EKF(3,3);
Q_EKF(4,4) = 100 .* Q_EKF(4,4);
Q_EKF(5,5) = 100 .* Q_EKF(5,5);
Q_EKF(6,6) = 1 .* Q_EKF(6,6);

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

    % Simulate noisy data from a randomly sampled initial state
    x_sim_init = mvnrnd(x_init, P_init)';
    [~, x_noise_mat, y_noise_mat] = simulateNoise(x_sim_init, u_nom, Q, R, dt, 1000);

    % Loop over the noisy measurements to get the estimates
    for i = 1:length(y_noise_mat(:,1))
        % Get the current measurement
        y_meas = y_noise_mat(i,:)';

        % Linearize the measurement function to get the Jacobian H
        [~, ~, H, ~] = linearize(xhat_curr, u_nom);

        % Kalman update
        [xhat_curr, P_curr, innovation, S] = EKF(xhat_curr, P_curr, u_nom, y_meas, dt, Q_EKF, R);

        % Store the result in a matrix of the estimated state
        xhat_mat(i+1, :) = xhat_curr';
        sigma_mat(i+1,:) = sqrt(diag(P_curr))';

        % Compute the NEES and NIS for this simulation
        % NEES (Normalized Estimation Error Squared)
        x_error = xhat_curr - x_noise_mat(i+1,:)';

        % Angle wrap the error
        x_error(3) = mod(x_error(3) + pi, 2*pi) - pi;
        x_error(6) = mod(x_error(6) + pi, 2*pi) - pi;

        % NEES
        nees_vals(sim_idx, i) = x_error' * inv(P_curr) * x_error;

        % NIS
        nis_vals(sim_idx, i) = innovation' * inv(S) * innovation;
    end
end

% Chi-Square Test Bounds
alpha = 0.05; % Significance level
r1_NEES = chi2inv(alpha/2, Nsim*Nstate) / Nsim;
r2_NEES = chi2inv(1-alpha/2, Nsim*Nstate) / Nsim;

r1_NIS = chi2inv(alpha/2, Nsim*Nmeas) / Nsim;
r2_NIS = chi2inv(1-alpha/2, Nsim*Nmeas) / Nsim;

% Average NEES and NIS across all runs
mean_nees = mean(nees_vals, 1); % Time-averaged NEES
mean_nis = mean(nis_vals, 1);   % Time-averaged NIS


% Plot the results
figure;
plot(time_tmt(2:end), mean_nees, 'ro', 'LineWidth', 1.5);
hold on;
yline(r1_NEES, 'r--', 'LineWidth', 1.2);
yline(r2_NEES, 'r--', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('NEES');
legend('Mean NEES', '\chi^2 Lower Bound', '\chi^2 Upper Bound');
title('NEES Chi-Square Test for EKF');
grid on;

figure
plot(time_tmt(2:end), mean_nis, 'bo', 'LineWidth', 1.5);
hold on;
yline(r1_NIS, 'b--', 'LineWidth', 1.2);
yline(r2_NIS, 'b--', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('NIS');
legend('Mean NIS', '\chi^2 Lower Bound', '\chi^2 Upper Bound');
title('NIS Chi-Square Test for EKF');
grid on;

%sgtitle('NEES and NIS Validation for Extended Kalman Filter');

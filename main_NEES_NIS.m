%% Clean
close all; clear; clc;

%% LFK
% Setup EOM
% Ode45 Constants
Data = load('cooplocalization_finalproj_KFdata.mat');
dt = 0.1;
tspan = [0 100];
t_nom = (dt:dt:tspan(2))';
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

% Nominal values
x_ugv = [10; 0; pi/2];
x_uav = [-60; 0; -pi/2];
u_ugv = [2; -pi/18];
u_uav = [12; pi/25];

x_nom = [x_ugv; x_uav];
u_nom = [u_ugv; u_uav];

% Simulate Nominal Nonlinear Trajectoy without Noise
w = zeros(6,1);
eomFunc = @(t, x)coopEOM(t, x, u_nom, w);
x_init = x_nom;
t_nom = (0:dt:tspan(2))';
[~, x_nom_mat] = ode45(eomFunc, t_nom, x_init, options);
u_nom_mat = ones(length(t_nom), 4) .* u_nom';

% Calculate the measurements from the sensor model
y_nom_mat = zeros(length(t_nom)-1, 5);
for i = 2:length(t_nom)
    y_nom_mat(i-1,:) = sensors(x_nom_mat(i,:))';
end

% Simulate Nonlinear Trajectory with Process Noise
Q_true = Data.Qtrue;
R_true = Data.Rtrue;

% Apply Linearized Kalman Filter
Q_tune = Q_true;
Q_tune(1,1) = Q_tune(1,1)*1000;
Q_tune(2,2) = Q_tune(2,2)*100;
Q_tune(3,3) = Q_tune(3,3)*100000;
Q_tune(4,4) = Q_tune(4,4)*100;
Q_tune(5,5) = Q_tune(5,5)*100;
Q_tune(6,6) = Q_tune(6,6)*10000;

% Monte Carlo Simulations for NEES and NIS
Nsim = 50; % Number of Monte Carlo simulations
Nstate = size(x_nom, 1);
Nmeas = size(y_nom_mat, 2);

% Preallocate for NEES and NIS results
nees_values = zeros(Nsim, length(t_nom)-1);
nis_values = zeros(Nsim, length(t_nom)-1);

for sim_idx = 1:Nsim
    % Simulate noisy trajectory
    [~, x_noisy, y_noisy] = simulateNoise(x_nom, u_nom, Q_true, R_true, dt, 1000);

    % Apply Linearized Kalman Filter
    [x_LKF, sigma, innovation_vec, S_vec] = LKF(x_nom_mat', u_nom_mat', y_nom_mat', y_noisy', u_nom_mat', Q_true, R_true, dt);

    % Compute NEES and NIS for this run
    for k = 1:length(t_nom)-1
        % State estimation error
        e_k = x_noisy(k, :)' - x_LKF(:, k);
        P_k = diag(sigma(:, k));

        % NEES (normalized state error)
        nees_values(sim_idx, k) = e_k' * (P_k \ e_k); % e_k' * inv(P_k) * e_k

        % NIS (normalized innovation error)
        innovation = innovation_vec(:,k);
        S = S_vec{k};
        nis_values(sim_idx, k) = innovation' * (S \ innovation); % innov' * inv(S_k) * innov
    end
end

% Chi-Square Test Bounds
alpha = 0.05; % Significance level
r1_NEES = chi2inv(alpha/2, Nsim*Nstate) / Nsim;
r2_NEES = chi2inv(1-alpha/2, Nsim*Nstate) / Nsim;

r1_NIS = chi2inv(alpha/2, Nsim*Nmeas) / Nsim;
r2_NIS = chi2inv(1-alpha/2, Nsim*Nmeas) / Nsim;

% Average NEES and NIS across all runs
mean_nees = mean(nees_values, 1); % Time-averaged NEES
mean_nis = mean(nis_values, 1);   % Time-averaged NIS


% Plotting
figure;
plot(t_nom(2:end), mean_nees, 'ro', 'LineWidth', 1.5);
hold on;
yline(r1_NEES, 'r--', 'LineWidth', 1.2);
yline(r2_NEES, 'r--', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('NEES');
legend('Mean NEES', '\chi^2 Lower Bound', '\chi^2 Upper Bound');
title('NEES Chi-Square Test for LKF');
grid on;

figure
plot(t_nom(2:end), mean_nis, 'bo', 'LineWidth', 1.5);
hold on;
yline(r1_NIS, 'b--', 'LineWidth', 1.2);
yline(r2_NIS, 'b--', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('NIS');
legend('Mean NIS', '\chi^2 Lower Bound', '\chi^2 Upper Bound');
title('NIS Chi-Square Test for LKF');
grid on;


%% EKF
clear; clc
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
Q = Q;
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
    [time_tmt, x_noise_mat, y_noise_mat] = simulateNoise(x_nom, u_nom, Q, R, dt, n_ind);
    xhat_curr = x_noise_mat(1,:);
    P_curr = P_init;

    % Loop over the noisy measurements to get the estimates
    for i = 1:length(y_noise_mat(:,1))
        % Get the current measurement
        y_meas = y_noise_mat(i,:)';

        % Linearize the measurement function to get the Jacobian H
        [~, ~, H, ~] = linearize(xhat_curr, u_nom);

        % Kalman update
        [xhat_curr, P_curr, innovation, S] = EKF(xhat_curr, P_curr, u_nom, y_meas, dt, Q, R);

        % Store the result in a matrix of the estimated state
        xhat_mat(i+1, :) = xhat_curr';
        sigma_mat(i+1,:) = sqrt(diag(P_curr))';

        % Compute the NEES and NIS for this simulation
        x_error = xhat_curr - x_noise_mat(i+1,:)';
        nees_vals(sim_idx, i) = x_error' * inv(P_curr) * x_error;
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

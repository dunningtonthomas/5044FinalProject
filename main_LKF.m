%% Main Function for Linear Kalman Filter
clc;
clear; 
close all; 

%% Setup EOM
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

%% Simulate Nominal Nonlinear Trajectoy without Noise
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


%% Simulate Nonlinear Trajectory with Process Noise
Q_true = Data.Qtrue;
R_true = Data.Rtrue;
[t_noise,x_noise_mat,y_noise_mat] = simulateNoise(x_nom,u_nom,Q_true,R_true,dt,1000);

%% Apply Linearized Kalman Filter
Q_tune = Q_true;
Q_tune(1,1) = Q_tune(1,1)*1000;
Q_tune(2,2) = Q_tune(2,2)*1000;
Q_tune(3,3) = Q_tune(3,3)*1000;
Q_tune(4,4) = Q_tune(4,4)*1000;
Q_tune(5,5) = Q_tune(5,5)*1000;
Q_tune(6,6) = Q_tune(6,6)*1000;

% % 1
% Q_tune(1,2) = Q_tune(1,2);
% Q_tune(2,1) = Q_tune(1,2);
% 
% Q_tune(1,3) = Q_tune(1,3)+1/100;
% Q_tune(3,1) = Q_tune(1,3);
% 
% Q_tune(1,4) = Q_tune(1,4)+1/100;
% Q_tune(4,1) = Q_tune(1,4);
% 
% Q_tune(1,5) = Q_tune(1,5)+1/100;
% Q_tune(5,1) = Q_tune(1,5);
% 
% Q_tune(1,6) = Q_tune(1,6);
% Q_tune(6,1) = Q_tune(1,6);
% 
% % 2
% Q_tune(2,3) = Q_tune(2,3)+1/100;
% Q_tune(3,2) = Q_tune(2,3);
% 
% Q_tune(2,4) = Q_tune(2,4)+1/1000;
% Q_tune(4,2) = Q_tune(2,4);
% 
% Q_tune(2,5) = Q_tune(2,5)+1/1000;
% Q_tune(5,2) = Q_tune(2,5);
% 
% Q_tune(2,6) = Q_tune(2,6)+1/1000;
% Q_tune(6,2) = Q_tune(2,6);
% % 3
% Q_tune(4,3) = Q_tune(4,3)+1/1000;
% Q_tune(3,4) = Q_tune(4,3);
% 
% Q_tune(5,3) = Q_tune(5,3)+1/1000;
% Q_tune(3,5) = Q_tune(5,3);
% 
% Q_tune(6,3) = Q_tune(6,3)+1/1000;
% Q_tune(3,6) = Q_tune(6,3);
% % 4
% Q_tune(4,5) = Q_tune(4,5);
% Q_tune(5,4) = Q_tune(4,5);
% 
% Q_tune(4,6) = Q_tune(4,6)+1/1000;
% Q_tune(6,4) = Q_tune(4,6);
% 
% Q_tune(5,6) = Q_tune(5,6)+1/1000;
% Q_tune(6,5) = Q_tune(5,6);
% 
% Q_tune = Q_tune*1000;

%% Monte Carlo Simulations for NEES and NIS
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
    [x_LKF, sigma] = LKF(x_nom_mat', u_nom_mat', y_nom_mat', y_noisy', u_nom_mat', Q_true, R_true, dt);

    % Compute NEES and NIS for this run
    for k = 1:length(t_nom)-1
        % State estimation error
        e_k = x_noisy(k, :)' - x_LKF(:, k);
        P_k = diag(sigma(:, k));

        % NEES (normalized state error)
        nees_values(sim_idx, k) = e_k' * (P_k \ e_k); % e_k' * inv(P_k) * e_k

        % Innovation and covariance
        innov = y_noisy(k, :)' - y_nom_mat(k, :)';

        % Use the linearize function to compute the observation Jacobian (C)
        [~, ~, H, ~] = linearize(x_LKF(:, k), u_nom); % Compute Jacobians
        S_k = R_true + H * P_k * H';

        % NIS (normalized innovation error)
        nis_values(sim_idx, k) = innov' * (S_k \ innov); % innov' * inv(S_k) * innov
    end
end

%% Chi-Square Test Bounds
alpha = 0.05; % Significance level
r1_NEES = chi2inv(alpha/2, Nstate) / Nsim;
r2_NEES = chi2inv(1-alpha/2, Nstate) / Nsim;

r1_NIS = chi2inv(alpha/2, Nmeas) / Nsim;
r2_NIS = chi2inv(1-alpha/2, Nmeas) / Nsim;

% Average NEES and NIS across all runs
mean_nees = mean(nees_values, 1); % Time-averaged NEES
mean_nis = mean(nis_values, 1);   % Time-averaged NIS


%% Plotting
% plotSim(t_noise, x_noise_mat, y_noise_mat, '--')
% plotSim(t_nom, x_nom_mat, y_nom_mat, '-')
% plotSim(t_noise, x_LKF', y_noise_mat, '-.')

x_error = x_LKF' - x_noise_mat;
% Angle wrapping
x_error(:,3) = mod(x_error(:,3) + pi, 2*pi) - pi;
x_error(:,6) = mod(x_error(:,6) + pi, 2*pi) - pi;
% State labels
state_labels = {'\xi_g Error [m]', '\eta_g Error [m]','\theta_g Error [rad]','\xi_a Error [m]','\eta_a Error [m]','\theta_a Error [rad]'};
% Plot the error for each state element of xs(k) vs time with ±2σ bounds
figure(1);
plot_num = 1;
for i = 1:6
    subplot(6, 1, plot_num);
    hold on;
    plot(t_noise, x_error(:,i), 'b', 'LineWidth', 1.5); 
    plot(t_noise(4:end), 2*sigma(i, 4:end), 'r--', 'LineWidth', 1.2);
    plot(t_noise(4:end), -2*sigma(i, 4:end), 'r--', 'LineWidth', 1.2);
    xlabel('Time [s]');
    ylabel(state_labels{plot_num});
    legend('Error', '\pm2\sigma', 'Location', 'best');
    title([state_labels{plot_num}, 'Error with \pm2\sigma Bounds']);
    grid on;
    plot_num =plot_num+1;
end
sgtitle('Aircrafts Position Estimator Error')

% Angle wrap before plotting
% x_LKF(3,:) = mod(x_LKF(3,:) + pi, 2*pi) - pi;
% x_LKF(6,:) = mod(x_LKF(6,:) + pi, 2*pi) - pi;
% 
% x_noise_mat(:,3) = mod(x_noise_mat(:,3) + pi, 2*pi) - pi;
% x_noise_mat(:,6) = mod(x_noise_mat(:,6) + pi, 2*pi) - pi;
% State labels
state_labels = {'\xi_g [m]', '\eta_g  [m]','\theta_g [rad]','\xi_a [m]','\eta_a  [m]','\theta_a [rad]'};
% Plot the error for each state element of xs(k) vs time with ±2σ bounds
figure(2);
plot_num = 1;
for i = 1:6
    subplot(6, 1, plot_num);
    hold on;
    plot(t_noise, x_LKF(i,:)', 'b', 'LineWidth', 1.5,'DisplayName','Estimated State'); 
    plot(t_noise, x_noise_mat(:,i)', 'r--', 'LineWidth', 1.2,'DisplayName','Nominal State');
    xlabel('Time [s]');
    ylabel(state_labels{plot_num});
    legend('Location', 'northeast');
    title(state_labels{plot_num});
    grid on;
    plot_num =plot_num+1;
end
sgtitle('Aircrafts State Over Time')


figure;
subplot(2, 1, 1);
plot(t_nom(2:end), mean_nees, 'b', 'LineWidth', 1.5);
hold on;
yline(r1_NEES, 'r--', 'LineWidth', 1.2);
yline(r2_NEES, 'r--', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('NEES');
legend('Mean NEES', '\chi^2 Lower Bound', '\chi^2 Upper Bound');
title('NEES Chi-Square Test');
grid on;

subplot(2, 1, 2);
plot(t_nom(2:end), mean_nis, 'b', 'LineWidth', 1.5);
hold on;
yline(r1_NIS, 'r--', 'LineWidth', 1.2);
yline(r2_NIS, 'r--', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('NIS');
legend('Mean NIS', '\chi^2 Lower Bound', '\chi^2 Upper Bound');
title('NIS Chi-Square Test');
grid on;

sgtitle('NEES and NIS Validation for Linearized Kalman Filter');
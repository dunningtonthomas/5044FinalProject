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

% Angle wrapping
x_nom_mat(:,3) = mod(x_nom_mat(:,3) + pi, 2*pi) - pi;
x_nom_mat(:,6) = mod(x_nom_mat(:,6) + pi, 2*pi) - pi;

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
Q_tune(1,1) = Q_tune(1,1)*100;
Q_tune(3,3) = Q_tune(3,3)*10;

[x_LKF,sigma]= LKF(x_nom_mat',u_nom_mat',y_nom_mat',y_noise_mat',u_nom_mat',Q_tune,R_true,dt);

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
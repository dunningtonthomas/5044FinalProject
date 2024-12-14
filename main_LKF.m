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
Q_tune(2,2) = Q_tune(2,2)*100;
Q_tune(3,3) = Q_tune(3,3)*100000;
Q_tune(4,4) = Q_tune(4,4)*100;
Q_tune(5,5) = Q_tune(5,5)*100;
Q_tune(6,6) = Q_tune(6,6)*10000;

% % % 1
% Q_tune(1,2) = Q_tune(1,2)+2.2/15;
% Q_tune(2,1) = Q_tune(1,2);
% 
% Q_tune(1,3) = Q_tune(1,3)-6;
% Q_tune(3,1) = Q_tune(1,3);
% 
% % % 2
% Q_tune(2,3) = Q_tune(2,3)+1.01;
% Q_tune(3,2) = Q_tune(2,3);
% 
% % 4
% Q_tune(4,5) = Q_tune(4,5)-1.3/20;
% Q_tune(5,4) = Q_tune(4,5);

% Q_tune(4,6) = Q_tune(4,6)-1/10;
% Q_tune(6,4) = Q_tune(4,6);
% 
% Q_tune(5,6) = Q_tune(5,6)-1/10;
% Q_tune(6,5) = Q_tune(5,6);
% 
% Q_tune = Q_tune*1000;
innovation = y_noise_mat-y_nom_mat;
% Angle wrap the innovation
% innovation(:,1) = mod(innovation(:,1) + pi, 2*pi) - pi;
% innovation(:,3) = mod(innovation(:,3) + pi, 2*pi) - pi;
% figure(3)
% hold on; 
% plot(t_noise(2:end),innovation);


[x_LKF,sigma,innovation_vec,S_vec]= LKF(x_nom_mat',u_nom_mat',y_nom_mat',y_noise_mat',u_nom_mat',Q_tune,R_true,dt);

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
    subplot(3, 2, plot_num);
    hold on;
    plot(t_noise, x_error(:,i), 'b', 'LineWidth', 1.5); 
    plot(t_noise(4:end), 2*sigma(i, 4:end), 'r--', 'LineWidth', 1.2);
    plot(t_noise(4:end), -2*sigma(i, 4:end), 'r--', 'LineWidth', 1.2);
    xlabel('Time [s]','FontSize',15);
    ylabel(state_labels{plot_num},'FontSize',15);
    legend('Error', '\pm2\sigma', 'Location', 'northwest','FontSize',20);
    title([state_labels{plot_num}, 'Error with \pm2\sigma Bounds'],'FontSize',20);
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
    subplot(3, 2, plot_num);
    hold on;
    plot(t_noise, x_LKF(i,:)', 'b', 'LineWidth', 1.5,'DisplayName','Estimated State'); 
    plot(t_noise, x_noise_mat(:,i)', 'r--', 'LineWidth', 1.2,'DisplayName','Nominal State');
    xlabel('Time [s]','FontSize',15);
    ylabel(state_labels{plot_num},'FontSize',15);
    legend('Location', 'northwest','FontSize',15);
    title(state_labels{plot_num},'FontSize',20);
    grid on;
    plot_num =plot_num+1;
end
sgtitle('Aircrafts State Over Time')
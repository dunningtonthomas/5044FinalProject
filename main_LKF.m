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
Q_tune(1,1) = Q_tune(1,1)*.1;
Q_tune(2,2) = Q_tune(2,2)*.1;
Q_tune(3,3) = Q_tune(3,3)*1000000;
Q_tune(4,4) = Q_tune(4,4)*0.1;
Q_tune(5,5) = Q_tune(5,5)*0.1;
Q_tune(6,6) = Q_tune(6,6)*0.00001;


Q_tune(1,2) = -6/1000000;
Q_tune(2,1) = Q_tune(1,2);

Q_tune(1,3) = 1/10000;
Q_tune(3,1) = Q_tune(1,3);

Q_tune(2,3) = 1/100000;
Q_tune(3,2) = Q_tune(2,3);

Q_tune(1,4) = 1/100000;
Q_tune(4,1) = Q_tune(1,4);

Q_tune(1,5) = 1/100000;
Q_tune(5,1) = Q_tune(1,5);


% Q_tune(4,5) = -1/10000;
% Q_tune(5,4) = Q_tune(4,5);



Q_tune = Q_tune.*200000;






% Q_tune = Q_true;
% Q_tune(1,1) = Q_tune(1,1)*1500;
% Q_tune(2,2) = Q_tune(2,2)*1000;
% Q_tune(3,3) = Q_tune(3,3)*8;
% Q_tune(4,4) = Q_tune(4,4)*1;
% Q_tune(5,5) = Q_tune(5,5)*1;
% Q_tune(6,6) = Q_tune(6,6)*.0001;
% Q_tune = Q_tune*200000;
% 
% Q_tune(1,2) = -.05;
% Q_tune(2,1) = Q_tune(1,2);
% 
% Q_tune(4,5) = -1/10;
% Q_tune(5,4) = Q_tune(4,5);

%% Run the LKF
[x_LKF,sigma,innovation_vec,S_vec]= LKF(x_nom_mat',u_nom_mat',y_nom_mat',y_noise_mat',u_nom_mat',Q_tune,R_true,dt);
%[x_LKF,sigma,innovation_vec,S_vec]= LKF(x_nom_mat',u_nom_mat',y_nom_mat',Data.ydata(:,2:end),u_nom_mat',Q_tune,R_true,dt);
%% Plotting
% plotSim(t_noise, x_noise_mat, y_noise_mat, '--')
% plotSim(t_nom, x_nom_mat, y_nom_mat, '-')
% plotSim(t_noise, x_LKF', y_noise_mat, '-.')

x_error = x_LKF' - x_noise_mat;
% Angle wrapping
% x_error(:,3) = mod(x_error(:,3) + pi, 2*pi) - pi;
% x_error(:,6) = mod(x_error(:,6) + pi, 2*pi) - pi;
% 
% x_LKF(3,:) = mod(x_LKF(3,:) + pi, 2*pi) - pi;
% x_LKF(6,:) = mod(x_LKF(6,:) + pi, 2*pi) - pi;
% 
% x_noise_mat(:,3) = mod(x_noise_mat(:,3) + pi, 2*pi) - pi;
% x_noise_mat(:,6) = mod(x_noise_mat(:,6) + pi, 2*pi) - pi;

%% Plot State Error with 2 Sigma Bound
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
    plot(t_noise, x_error(:,i)', 'b', 'LineWidth', 1.5, 'DisplayName', 'Error'); 
    plot(t_noise(4:end), 2*sigma(i,4:end), 'r--', 'LineWidth', 1.2, 'DisplayName', '+2\sigma'); 
    plot(t_noise(4:end), -2*sigma(i,4:end), 'r--', 'LineWidth', 1.2, 'DisplayName', '-2\sigma');
    % Add labels and grid
    xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel(state_labels{plot_num}, 'Interpreter', 'latex', 'FontSize', 15);
    % Add legend for the top right
    if(plot_num == 4)
        legend('Estimate', '$2 \sigma$ bounds', 'interpreter', 'latex');
    end
    % Add title
    title([state_labels{plot_num}, ' with $\pm 2 \sigma$ Bounds'], 'Interpreter', 'latex', 'FontSize', 15);
    % Increment plot number
    plot_num = plot_num + 1;
end

% Add a global title
sgtitle('State Estimator Errors with $\pm 2 \sigma$  Bounds', 'Interpreter', 'latex', 'FontSize', 20);
% Global adjustments
set(gcf, 'Position', [100, 100, 1200, 800]); % Adjust figure size


%% Plot Estimated State Vs Ground Truth w/ Process Noise
% State labels
state_labels = {'$\xi_g$ [m]', '$\eta_g$ [m]', '$\theta_g$ [rad]', ...
                '$\xi_a$ [m]', '$\eta_a$ [m]', '$\theta_a$ [rad]'};
% Create figure
figure('Name', 'Ground Truth States vs Estimated States', 'NumberTitle', 'off', 'Color', 'w');
plot_num = 1;
% Loop to create 3x2 subplot
plot_locations = [1, 3, 5, 2, 4, 6]; % Define plot order for better arrangement
for i = 1:6
    subplot(3, 2, plot_locations(plot_num)); % 3 rows, 2 columns
    hold on;
    % Plot estimated state and noisy ground truth state
    plot(t_noise, x_LKF(i, :)', 'b', 'LineWidth', 1.5, 'DisplayName', 'Estimated State'); 
    plot(t_noise, x_noise_mat(:, i)', 'r--', 'LineWidth', 1.2, 'DisplayName', 'Noisy Ground Truth State');
    % Add labels and grid
    xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel(state_labels{plot_num}, 'Interpreter', 'latex', 'FontSize', 15);
    % Add title
    title([state_labels{plot_num}], 'Interpreter', 'latex', 'FontSize', 15);
    % Add legend for the top-right plot
    if(plot_num == 4)
        legend('Estimated State', 'Noisy Ground Truth', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 12);
    end
    grid on;
    % Increment plot number
    plot_num = plot_num + 1;
end
% Add a global title
sgtitle('Ground Truth States vs Estimated States Over Time', 'Interpreter', 'latex', 'FontSize', 20);
% Global adjustments
set(gcf, 'Position', [100, 100, 1200, 800]); % Adjust figure size




%% Plot Estimated States with 2 Sigma Bounds
% % Create figure
% figure('Name', 'LKF State Estimates', 'NumberTitle', 'off', 'Color', 'w');
% plot_num = 1;
% % Loop to create 3x2 subplot
% plot_locations = [1, 3, 5, 2, 4, 6];
% for i = 1:6
%     subplot(3, 2, plot_locations(plot_num)); % 3 rows, 2 columns
%     hold on;
%     % Plot error and bounds
%     plot(t_noise, x_LKF(i,:)', 'b', 'LineWidth', 1.5, 'DisplayName', 'Estimate'); 
%     plot(t_noise,  x_LKF(i,:)+2*sigma(i,:), 'r--', 'LineWidth', 1.2, 'DisplayName', '+2\sigma'); 
%     plot(t_noise,  x_LKF(i,:)-2*sigma(i,:), 'r--', 'LineWidth', 1.2, 'DisplayName', '-2\sigma');
%     % Add labels and grid
%     xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 15);
%     ylabel(state_labels{plot_num}, 'Interpreter', 'latex', 'FontSize', 15);
%     % Add legend for the top right
%     if(plot_num == 4)
%         legend('Estimate', '$2 \sigma$ bounds', 'interpreter', 'latex');
%     end
%     % Add title
%     title([state_labels{plot_num}, ' with $\pm 2 \sigma$ Bounds'], 'Interpreter', 'latex', 'FontSize', 15);
%     % Increment plot number
%     plot_num = plot_num + 1;
% end
% 
% % Add a global title
% sgtitle('LKF State Estimates with $\pm 2 \sigma$  Bounds', 'Interpreter', 'latex', 'FontSize', 20);
% % Global adjustments
% set(gcf, 'Position', [100, 100, 1200, 800]); % Adjust figure size

%% Plot Innovation
% State labels
state_labels = {'y1 [rad]', 'y2 [m]', 'y3 [rad]', ...
                '$\xi_a$ [m]', '$\eta_a$ [m]'};
% Create figure
figure('Name', 'Innovation', 'NumberTitle', 'off', 'Color', 'w');
plot_num = 1;
% Loop to create 3x2 subplot
plot_locations = [1,2,3,4,5]; % Define plot order for better arrangement
for i = 1:5
    subplot(5, 1, plot_locations(plot_num)); % 3 rows, 2 columns
    hold on;
    % Plot estimated state and noisy ground truth state
    plot(t_noise(2:end), innovation_vec(i, :)', 'b', 'LineWidth', 1.5, 'DisplayName', 'Innovation'); 
    % Add labels and grid
    xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel(state_labels{plot_num}, 'Interpreter', 'latex', 'FontSize', 15);
    % Add title
    title([state_labels{plot_num}], 'Interpreter', 'latex', 'FontSize', 15);
    % Add legend for the top-right plot
    if(plot_num == 4)
        legend('Innovation', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 12);
    end
    grid on;
    % Increment plot number
    plot_num = plot_num + 1;
end
% Add a global title
sgtitle('Innovation Over Time', 'Interpreter', 'latex', 'FontSize', 20);
% Global adjustments
set(gcf, 'Position', [100, 100, 1200, 800]); % Adjust figure size

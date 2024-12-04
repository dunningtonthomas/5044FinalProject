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
n_ind = length(t_nom)-1;
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
y_nom_mat = zeros(length(t_nom), 5);
for i = 1:length(t_nom)
    y_nom_mat(i,:) = sensors(x_nom_mat(i,:))';
end


%% Simulate Nonlinear Trajectory with Process Noise
rng(100);
Q_true = Data.Qtrue;
Sv_process = chol(Q_true,'lower');
x_noise_mat = x_nom';
t_noise = 0;
for i=1:n_ind+1
    % Define time step
    ZOH_time = dt*[i-1 i];
    % Generate AWGN for the time step
    q = randn(6,1);
    w =  Sv_process*q;
    % Set Initial Conditions
    x_init = x_noise_mat(end,:);
    % Function Handle for the nominal trajectory using the full nonlinear equations with noise
    eomFunc = @(t, x)coopEOM(t, x, u_nom, w);
    % Simulate over time step
    [TOUT, x_noise_step] = ode45(eomFunc, ZOH_time, x_init, options);
    x_noise_mat = [x_noise_mat;x_noise_step(end,:)];
    t_noise = [t_noise;TOUT(end)];
end

% Angle wrapping
x_noise_mat(:,3) = mod(x_noise_mat(:,3) + pi, 2*pi) - pi;
x_noise_mat(:,6) = mod(x_noise_mat(:,6) + pi, 2*pi) - pi;

% Calculate the measurements from the sensor model
y_noise_mat = zeros(length(t_noise), 5);
for i = 1:length(t_noise)
    y_noise_mat(i,:) = sensors(x_noise_mat(i,:))';
end
% Add Measurement Noise 
R_true = Data.Rtrue;
Sv_measurement = chol(R_true,'lower');
q = randn(5,length(t_noise));
y_noise_mat =  y_noise_mat + (Sv_measurement*q)';

%% Apply Linearized Kalman Filter
[x_LKF,sigma]= LKF(x_nom_mat',u_nom_mat',y_nom_mat',y_noise_mat',u_nom_mat',Q_true,R_true,dt);


%% Plotting
% plotSim(t_noise, x_noise_mat, y_noise_mat, '--')
plotSim(t_nom, x_nom_mat, y_nom_mat, '-')
plotSim(t_noise(2:end), x_LKF(:,2:end)', y_noise_mat(2:end,:), '-.')
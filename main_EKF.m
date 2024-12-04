%% Main Function for EKF Analysis
close all; clear; clc;

%% EKF
Data = load('cooplocalization_finalproj_KFdata.mat');

% Test the EKF update
% Nominal values
x_ugv = [10; 0; pi/2];
x_uav = [-60; 0; -pi/2];
u_ugv = [2; -pi/18];
u_uav = [12; pi/25];

x_nom = [x_ugv; x_uav];
u_nom = [u_ugv; u_uav];

xhat_prev = x_nom;
u = u_nom;
y = sensors(xhat_prev);
dt = 0.1;
Q = Data.Qtrue;
R = Data.Rtrue;
P_prev = [10, 0, 0, 0, 0, 0;
            0, 10, 0, 0, 0, 0;
            0, 0, pi, 0, 0, 0;
            0, 0, 0, 10, 0, 0;
            0, 0, 0, 0, 10, 0;
            0, 0, 0, 0, 0, pi];

[xhat_meas, P_meas] = EKF(xhat_prev, P_prev, u, y, dt, Q, R);







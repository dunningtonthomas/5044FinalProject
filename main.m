%% MAIN FUNCTION
close all; clear; clc;

%% Simulate EOM
% Nominal values
x_ugv = [10; 0; pi/2];
x_uav = [-60; 0; -pi/2];
u_ugv = [2; -pi/18];
u_uav = [12; pi/25];

x_nom = [x_ugv; x_uav];
u_nom = [u_ugv; u_uav];

% Simulate EOM
eomFunc = @(t, x)coopEOM(t, x, u_nom);
x_init = x_nom;
tspan = [0 50];
[TOUT, YOUT] = ode45(eomFunc, tspan, x_init);


% Plot
plotSim(TOUT, YOUT)


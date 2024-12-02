%% MAIN FUNCTION
close all; clear; clc;

%% Simulate EOM
% Constants
dt = 0.1;
tspan = [0 50];

% Nominal values
x_ugv = [10; 0; pi/2];
x_uav = [-60; 0; -pi/2];
u_ugv = [2; -pi/18];
u_uav = [12; pi/25];

x_nom = [x_ugv; x_uav];
u_nom = [u_ugv; u_uav];

% Linearize about the nominal point
[A, B, C, D] = linearize(x_nom, u_nom);

% Find the DT model from the linearized CT model
[F, G, H, M] = discretize(A, B, C, D, dt);

% Find the nominal trajectory using the full nonlinear equations
eomFunc = @(t, x)coopEOM(t, x, u_nom);
x_init = x_nom;
times = (dt:dt:tspan(2))';
[~, x_nom_mat] = ode45(eomFunc, times, x_init);
u_nom_mat = ones(length(times), 4) .* u_nom';


% Simulate the discrete model with an initial perturbation
dx0 = [0.1; 0.1; 0.1; 0.1; 0.1; 0.1];
dx0 = [0; 0; 0; 0; 0; 0];
[XOUT_DT, YOUT_DT] = simulateDT(x_nom_mat, u_nom_mat, dx0, times);
plotSim(times, XOUT_DT)

% Plot
%plotSim(TOUT, YOUT)



%% MAIN FUNCTION
close all; clear; clc;

%% Simulate EOM
% Ode45 Constants
dt = 0.1;
tspan = [0 100];
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

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

% Find the nominal trajectory using the full nonlinear equations with no
% noise
w = zeros(6,1);
eomFunc = @(t, x)coopEOM(t, x, u_nom, w);
x_init = x_nom;
times = (dt:dt:tspan(2))';
[~, x_nom_mat] = ode45(eomFunc, times, x_init, options);
u_nom_mat = ones(length(times), 4) .* u_nom';


% Simulate the discrete model with an initial perturbation
dx0 = [0; 1; 0; 0; 0; 0.1];
[XOUT_DT, YOUT_DT] = simulateDT(x_nom_mat, u_nom_mat, dx0, times);

% Simulate full nonlinear EOM
x_init = x_nom + dx0;
[TOUT_NL, XOUT_NL] = ode45(eomFunc, tspan, x_init, options);

% Angle wrapping
XOUT_NL(:,3) = mod(XOUT_NL(:,3) + pi, 2*pi) - pi;
XOUT_NL(:,6) = mod(XOUT_NL(:,6) + pi, 2*pi) - pi;

% Calculate the measurements from the sensor model
YOUT_NL = zeros(length(TOUT_NL), 5);
for i = 1:length(TOUT_NL)
    YOUT_NL(i,:) = sensors(XOUT_NL(i,:))';
end


% Plot
%plotSim(times, XOUT_DT, '-')
plotSim(TOUT_NL, XOUT_NL, YOUT_NL, '-')



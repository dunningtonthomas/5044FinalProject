function [XOUT, YOUT] = simulateDT(x_nom_mat, u_nom_mat, dx0, times)
%SIMULATEDT This function simulates the DT linear model and returns the
%time, states, and measurements
%
% Inputs: 
%   x_nom_mat -> nominal states [length(t) x 6]
%   u_nom_mat -> nominal trajectory inputs [length(t)x4]
%   dx0 -> initial state perturbation
%   times -> vector of discrete time values to propagate the dynamics at
%
% Outputs:
%   XOUT -> Full discrete time state output
%   YOUT -> Full discrete time measurements output
%
% Author: Thomas Dunnington
% Modified: 12/2/2024

% Discrete time step
dt = times(2) - times(1);

% Initial state and measurement
XOUT = zeros(length(times), 6);
xk = zeros(length(times), 6);
YOUT = zeros(length(times), 5);
yk = zeros(length(times), 5);

% Initial state
xk(1, :) = dx0;
XOUT(1, :) = x_nom_mat(1,:) + xk(1, :);

% Simulate state perturbations over time
for i = 1:length(times)-1
    % Get current time and nominal state/control
    time = times(i);
    x_nom = x_nom_mat(i,:)';
    u_nom = u_nom_mat(i,:)';

    % Linearize around the nominal point
    [A, B, C, D] = linearize(x_nom, u_nom);

    % Discrete time conversion
    [F, G, H, M] = discretize(A, B, C, D, dt);

    % Simulate DT dynamics
    xk(i+1, :) = (F*xk(i, :)')';
    yk(i, :) = (H*xk(i, :)')';

    % Add the perturbation to the nominal trajectory
    XOUT(i+1, :) = xk(i+1, :) + x_nom_mat(i+1,:);
end

% Get the last measurement
x_nom = x_nom_mat(end,:)';
u_nom = u_nom_mat(end,:)';

% Linearize around the nominal point
[A, B, C, D] = linearize(x_nom, u_nom);

% Discrete time conversion
[F, G, H, M] = discretize(A, B, C, D, dt);

% Simulate DT dynamics
YOUT(end, :) = (H*xk(end,:)')';

end


%% MAIN FUNCTION PART 1
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
TOUT_DT = (0:dt:tspan(2))';
[~, x_nom_mat] = ode45(eomFunc, TOUT_DT, x_init, options);
u_nom_mat = ones(length(TOUT_DT), 4) .* u_nom';


% Simulate the discrete model with an initial perturbation
dx0 = [0; 1; 0; 0; 0; 0.1];
[XOUT_DT, YOUT_DT] = simulateDT(x_nom_mat, u_nom_mat, dx0, TOUT_DT);

% Simulate full nonlinear EOM
x_init = x_nom + dx0;
[TOUT_NL, XOUT_NL] = ode45(eomFunc, TOUT_DT, x_init, options);

% Calculate the measurements for the NL simulation
YOUT_NL = zeros(length(TOUT_NL)-1, 5);
for i = 2:length(TOUT_NL)
    YOUT_NL(i-1,:) = sensors(XOUT_NL(i,:))';
end

% Plot
%plotSim(times, XOUT_DT, '-')
plotSim(TOUT_NL, XOUT_NL, YOUT_NL, '-')
%plotSim(TOUT_DT, XOUT_DT, YOUT_DT, '--')

% Plot the perturbations of the linear approximation
perturbationPlot(TOUT_DT, XOUT_DT, x_nom_mat, '-');


% Calculate and plot the difference between the simulated states and the
% measurements
xdiff = XOUT_DT - XOUT_NL;
ydiff = YOUT_DT - YOUT_NL;


% Plot the difference between the linear and nonlinear states
figure();
sgtitle('Difference Between Nonlinear and Linear Approximation', 'FontSize', 16, 'FontWeight', 'bold'); % Enhanced title appearance

% Define line properties for reusability
line_width = 1.5;
linespec = '-';
ugv_line = {'LineWidth', line_width, 'Color', 'r', 'LineStyle', linespec};
uav_line = {'LineWidth', line_width, 'Color', 'r', 'LineStyle', linespec};

% UGV East Position
subplot(3,2,1);
plot(TOUT_DT, xdiff(:,1), ugv_line{:});
ylabel('$\Delta\xi_{g}$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
ylim([-0.25 0.25])
title('UGV States', 'FontSize', 12);
grid on;

% UAV East Position
subplot(3,2,2);
plot(TOUT_DT, xdiff(:,4), uav_line{:});
ylabel('$\Delta\xi_{a}$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
title('UAV States', 'FontSize', 12);
grid on;

% UGV North Position
subplot(3,2,3);
plot(TOUT_DT, xdiff(:,2), ugv_line{:});
ylim([-0.25 0.25])
ylabel('$\Delta\eta_{g}$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
grid on;

% UAV North Position
subplot(3,2,4);
plot(TOUT_DT, xdiff(:,5), uav_line{:});
ylabel('$\Delta\eta_{a}$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
grid on;

% UGV Heading Angle
subplot(3,2,5);
plot(TOUT_DT, xdiff(:,3), ugv_line{:});
ylim([-0.25 0.25])
ylabel('$\Delta\theta_{g}$ (rad)', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
grid on;

% UAV Heading Angle
subplot(3,2,6);
plot(TOUT_DT, xdiff(:,6), uav_line{:});
ylim([-0.25 0.25])
ylabel('$\Delta\theta_{a}$ (rad)', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
grid on;

% Global adjustments
set(gcf, 'Position', [100, 100, 1200, 800]); % Adjust figure size


% Difference in the measurements
figure();
sgtitle('Difference Between Linear and Nonlinear Measurements', 'FontSize', 14, 'FontWeight', 'bold');

% Define y-axis labels for each subplot
ylabels = {
    '$\Delta\gamma_{ag}$ (rad)', ...
    '$\Delta\rho_{ga}$ (m)', ...
    '$\Delta\gamma_{ga}$ (rad)', ...
    '$\Delta\xi_{a}$ (m)', ...
    '$\Delta\eta_{a}$ (m)'
};

% Loop through subplots
for i = 1:5
    subplot(5, 1, i);
    plot(TOUT_DT(2:end), ydiff(:, i), 'LineWidth', 1.5, 'Color', 'r', 'LineStyle', linespec); % Set line properties
    ylabel(ylabels{i}, 'Interpreter', 'latex', 'FontSize', 12); % Y-axis label
    grid on; % Add grid
    hold on;
    ytickformat('%.2f'); % Simplify y-axis tick format
    if i < 5
        set(gca, 'XTickLabel', []); % Remove x-tick labels for all but the last subplot
    else
        xlabel('Time (s)', 'FontSize', 12); % Add x-axis label only on the last subplot
    end
end

% Adjust subplot spacing for better clarity
set(gcf, 'Position', [100, 100, 800, 600]); % Resize figure for better proportions




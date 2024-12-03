function perturbationPlot(TOUT, XOUT, XNOM, linespec)
%PERTURBATIONPLOT Plot the perturbations of each state given a nominal
%trajectory and the actual trajectory
%states over time
% 
% Inputs:
%   TOUT -> time vector
%   XOUT -> length(TOUT) by 6 state matrix, each column is a state
%   XNOM -> nominal trajectory matrix
%   linespec -> linestyle of the plot
%
% Output:
%   plot 1 -> 6x1 subplot of the state perturbations
%
% Author: Thomas Dunnington
% Modified: 12/3/2024

% Calculate the perturbations
XPERT = XOUT - XNOM;

% figure();
% sgtitle('Linearized State Perturbations', 'FontSize', 14, 'FontWeight', 'bold');
% subplot(6,1,1)
% plot(TOUT, XPERT(:,1), 'linewidth', 1.5, 'color', [0, 0.4470, 0.7410]);
% ylabel('$\delta\xi_{g}$ (m)', 'Interpreter', 'latex')
% grid on;
% 
% 
% subplot(6,1,2)
% plot(TOUT, XPERT(:,2), 'linewidth', 1.5, 'color', [0, 0.4470, 0.7410]);
% ylabel('$\delta\eta_{g}$ (m)', 'Interpreter', 'latex')
% grid on;
% ylim([0 2])
% 
% 
% subplot(6,1,3)
% plot(TOUT, XPERT(:,3), 'linewidth', 1.5, 'color', [0, 0.4470, 0.7410]);
% ylabel('$\delta\theta_{g}$ (rad)', 'Interpreter', 'latex')
% grid on;
% 
% subplot(6,1,4)
% plot(TOUT, XPERT(:,4), 'linewidth', 1.5, 'color', [0, 0.4470, 0.7410]);
% ylabel('$\delta\xi_{a}$ (m)', 'Interpreter', 'latex')
% grid on;
% 
% subplot(6,1,5)
% plot(TOUT, XPERT(:,5), 'linewidth', 1.5, 'color', [0, 0.4470, 0.7410]);
% ylabel('$\delta\eta_{a}$ (m)', 'Interpreter', 'latex')
% grid on;
% 
% subplot(6,1,6)
% plot(TOUT, XPERT(:,6), 'linewidth', 1.5, 'color', [0, 0.4470, 0.7410]);
% ylabel('$\delta\theta_{a}$ (rad)', 'Interpreter', 'latex')
% grid on;
% ylim([-1 1])
% xlabel('Time (s)')
% 
% set(gcf, 'Position', [100, 100, 800, 600]); % Resize figure for better proportions



figure();
sgtitle('Linearized State Perturbations', 'FontSize', 16, 'FontWeight', 'bold');

lineColor = [0, 0.4470, 0.7410]; % Consistent blue color for all plots

subplot(6,1,1)
plot(TOUT, XPERT(:,1), 'LineWidth', 1.5, 'Color', lineColor);
ylabel('$\delta\xi_{g}$ (m)', 'Interpreter', 'latex', 'FontSize', 12)
grid on;

subplot(6,1,2)
plot(TOUT, XPERT(:,2), 'LineWidth', 1.5, 'Color', lineColor);
ylabel('$\delta\eta_{g}$ (m)', 'Interpreter', 'latex', 'FontSize', 12)
grid on;
ylim([0 2]);

subplot(6,1,3)
plot(TOUT, XPERT(:,3), 'LineWidth', 1.5, 'Color', lineColor);
ylabel('$\delta\theta_{g}$ (rad)', 'Interpreter', 'latex', 'FontSize', 12)
grid on;

subplot(6,1,4)
plot(TOUT, XPERT(:,4), 'LineWidth', 1.5, 'Color', lineColor);
ylabel('$\delta\xi_{a}$ (m)', 'Interpreter', 'latex', 'FontSize', 12)
grid on;

subplot(6,1,5)
plot(TOUT, XPERT(:,5), 'LineWidth', 1.5, 'Color', lineColor);
ylabel('$\delta\eta_{a}$ (m)', 'Interpreter', 'latex', 'FontSize', 12)
grid on;

subplot(6,1,6)
plot(TOUT, XPERT(:,6), 'LineWidth', 1.5, 'Color', lineColor);
ylabel('$\delta\theta_{a}$ (rad)', 'Interpreter', 'latex', 'FontSize', 12)
grid on;
ylim([-1 1]);
xlabel('Time (s)', 'FontSize', 12)

% Adjust spacing and figure size
set(gcf, 'Position', [100, 100, 1000, 800]); % Resize figure for better proportions



end


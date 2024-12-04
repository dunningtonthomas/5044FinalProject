function plotSim(TOUT, XOUT, YOUT, linespec)
%PLOTSIM Plots the 2D positions of the UGV and UAV and a 3x2 subplot of the
%states over time
% 
% Inputs:
%   TOUT -> time vector
%   XOUT -> length(TOUT) by 6 state matrix, each column is a state
%   YOUT -> length(TOUT)-1 by 5 measurement matrix
%   linespec -> linestyle of the plot
%
% Output:
%   plot 1 -> 2D positions of UAV and UGV
%   plot 2 -> 3x2 subplot of each state over time
%   plot 3 -> 5x1 subplot of the measurements with no measurement noise
%
% Author: Thomas Dunnington
% Modified: 12/3/2024

% Angle wrapping
XOUT(:,3) = mod(XOUT(:,3) + pi, 2*pi) - pi;
XOUT(:,6) = mod(XOUT(:,6) + pi, 2*pi) - pi;

% Define colors for better visualization
ugv_color = [0, 0.4470, 0.7410]; % UGV color
uav_color = [0.8500, 0.3250, 0.0980]; % UAV color

% Enhanced comparison plot
figure(1);
plot(XOUT(:,1), XOUT(:,2), 'LineWidth', 1.5, 'Color', ugv_color, 'DisplayName', 'UGV', 'linestyle', linespec)
hold on;
plot(XOUT(:,4), XOUT(:,5), 'LineWidth', 1.5, 'Color', uav_color, 'DisplayName', 'UAV', 'linestyle', linespec)
xlabel('East Position');
ylabel('North Position');
title('UGV vs UAV Position');
legend('Location', 'best');
grid on;
axis on;

% Enhanced subplots
figure(2);
sgtitle('Simulated States')
subplot(3,2,1)
plot(TOUT, XOUT(:,1), 'LineWidth', 1.5, 'Color', ugv_color, 'linestyle', linespec)
ylabel('East Position');
title('UGV East Position');
grid on;
hold on;

subplot(3,2,2)
plot(TOUT, XOUT(:,4), 'LineWidth', 1.5, 'Color', uav_color, 'linestyle', linespec)
ylabel('East Position');
title('UAV East Position');
grid on;
hold on;

subplot(3,2,3)
plot(TOUT, XOUT(:,2), 'LineWidth', 1.5, 'Color', ugv_color, 'linestyle', linespec)
ylabel('North Position');
title('UGV North Position');
grid on;
hold on;

subplot(3,2,4)
plot(TOUT, XOUT(:,5), 'LineWidth', 1.5, 'Color', uav_color, 'linestyle', linespec)
ylabel('North Position');
title('UAV North Position');
grid on;
hold on;

subplot(3,2,5)
plot(TOUT, XOUT(:,3), 'LineWidth', 1.5, 'Color', ugv_color, 'linestyle', linespec)
ylabel('Heading Angle');
title('UGV Heading Angle');
xlabel('Time');
grid on;
hold on;

subplot(3,2,6)
plot(TOUT, XOUT(:,6), 'LineWidth', 1.5, 'Color', uav_color, 'linestyle', linespec)
ylabel('Heading Angle');
title('UAV Heading Angle');
xlabel('Time');
grid on;
hold on;

% Global figure adjustments
set(gcf, 'Position', [100, 100, 1200, 800]); % Adjust figure size


% Measurement vector plot
figure(3);
sgtitle('Measurements', 'FontSize', 14, 'FontWeight', 'bold');

% Define y-axis labels for each subplot
ylabels = {
    '$\gamma_{ag}$ (rad)', ...
    '$\rho_{ga}$ (m)', ...
    '$\gamma_{ga}$ (rad)', ...
    '$\xi_{a}$ (m)', ...
    '$\eta_{a}$ (m)'
};

% Loop through subplots
for i = 1:5
    subplot(5, 1, i);
    plot(TOUT(2:end), YOUT(:, i), 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', linespec); % Set line properties
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



end


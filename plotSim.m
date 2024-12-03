function plotSim(TOUT, XOUT, YOUT, linespec)
%PLOTSIM Plots the 2D positions of the UGV and UAV and a 3x2 subplot of the
%states over time
% 
% Inputs:
%   TOUT -> time vector
%   XOUT -> length(TOUT) by 6 state matrix, each column is a state
%   YOUT -> length(TOUT) by 5 measurement matrix
%   linespec -> linestyle of the plot
%
% Output:
%   plot 1 -> 2D positions of UAV and UGV
%   plot 2 -> 3x2 subplot of each state over time
%   plot 3 -> 5x1 subplot of the measurements with no measurement noise

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


% Measurement vector
figure(3);
sgtitle('Measurement')
subplot(5,1,1)
plot(TOUT, YOUT(:,1), 'linewidth', 1.5, 'Color', 'k', 'linestyle', linespec)
ylabel('$\gamma_{ag}$ (rad)', 'Interpreter', 'latex')

subplot(5,1,2)
plot(TOUT, YOUT(:,2), 'linewidth', 1.5, 'Color', 'k', 'linestyle', linespec)
ylabel('$\rho_{ga}$ (m)', 'Interpreter', 'latex')

subplot(5,1,3)
plot(TOUT, YOUT(:,3), 'linewidth', 1.5, 'Color', 'k', 'linestyle', linespec)
ylabel('$\gamma_{ga}$ (rad)', 'Interpreter', 'latex')

subplot(5,1,4)
plot(TOUT, YOUT(:,4), 'linewidth', 1.5, 'Color', 'k', 'linestyle', linespec)
ylabel('$\xi_{a}$ (m)', 'Interpreter', 'latex')

subplot(5,1,5)
plot(TOUT, YOUT(:,5), 'linewidth', 1.5, 'Color', 'k', 'linestyle', linespec)
ylabel('$\eta_{a}$ (m)', 'Interpreter', 'latex')
xlabel('Time (s)')


end


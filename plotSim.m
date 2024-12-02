function plotSim(TOUT, YOUT, linespec)
%PLOTSIM Plots the 2D positions of the UGV and UAV and a 3x2 subplot of the
%states over time
% 
% Inputs:
%   TOUT -> time vector
%   YOUT -> length(TOUT) by 6 state matrix, each column is a state over
%   linespec -> linestyle of the plot
%
% Output:
%   plot 1 -> 2D positions of UAV and UGV
%   plot 2 -> 3x2 subplot of each state over time

% Define colors for better visualization
ugv_color = [0, 0.4470, 0.7410]; % UGV color
uav_color = [0.8500, 0.3250, 0.0980]; % UAV color

% Enhanced comparison plot
figure(1);
plot(YOUT(:,1), YOUT(:,2), 'LineWidth', 1.5, 'Color', ugv_color, 'DisplayName', 'UGV', 'linestyle', linespec)
hold on;
plot(YOUT(:,4), YOUT(:,5), 'LineWidth', 1.5, 'Color', uav_color, 'DisplayName', 'UAV', 'linestyle', linespec)
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
plot(TOUT, YOUT(:,1), 'LineWidth', 1.5, 'Color', ugv_color, 'linestyle', linespec)
ylabel('East Position');
title('UGV East Position');
grid on;
hold on;

subplot(3,2,2)
plot(TOUT, YOUT(:,4), 'LineWidth', 1.5, 'Color', uav_color, 'linestyle', linespec)
ylabel('East Position');
title('UAV East Position');
grid on;
hold on;

subplot(3,2,3)
plot(TOUT, YOUT(:,2), 'LineWidth', 1.5, 'Color', ugv_color, 'linestyle', linespec)
ylabel('North Position');
title('UGV North Position');
grid on;
hold on;

subplot(3,2,4)
plot(TOUT, YOUT(:,5), 'LineWidth', 1.5, 'Color', uav_color, 'linestyle', linespec)
ylabel('North Position');
title('UAV North Position');
grid on;
hold on;

subplot(3,2,5)
plot(TOUT, YOUT(:,3), 'LineWidth', 1.5, 'Color', ugv_color, 'linestyle', linespec)
ylabel('Heading Angle');
title('UGV Heading Angle');
xlabel('Time');
grid on;
hold on;

subplot(3,2,6)
plot(TOUT, YOUT(:,6), 'LineWidth', 1.5, 'Color', uav_color, 'linestyle', linespec)
ylabel('Heading Angle');
title('UAV Heading Angle');
xlabel('Time');
grid on;
hold on;

% Global figure adjustments
set(gcf, 'Position', [100, 100, 1200, 800]); % Adjust figure size


end


function xdot = coopEOM(t, x, u)
%COOPEOM Full nonlinear equations of motion for the UGV and UAV with no
%process noise
%
% Inputs:
%   t -> time
%   x -> Full state (xi_g, eta_g, theta_g, xi_a, eta_a, theta_a)
%   u -> control inputs (v_g, phi_g, v_a, w_a)
%
% Outputs:
%   xdot -> time rate of change of the state
%
% Author: Thomas Dunnington
% Modified: 12/2/2024

% Separation length
L = 0.5;

% UGV EOM
v_g = u(1); theta_g = x(3); phi_g = u(2);
xi_g_dot = v_g*cos(theta_g);
eta_g_dot = v_g*sin(theta_g);
theta_g_dot = v_g/L * tan(phi_g);
x_ugv_dot = [xi_g_dot; eta_g_dot; theta_g_dot];

% UAV EOM
v_a = u(3); w_a = u(4); theta_a = x(6);
xi_a_dot = v_a*cos(theta_a);
eta_a_dot = v_a*sin(theta_a);
theta_a_dot = w_a;
x_uav_dot = [xi_a_dot; eta_a_dot; theta_a_dot];

% Full state ROC
xdot = [x_ugv_dot; x_uav_dot];

end


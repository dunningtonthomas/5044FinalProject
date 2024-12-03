function [A, B, C, D] = linearize(x_nom, u_nom)
%LINEARIZE This function linearizes the nonlinear equations by calculating
%the jacobians around the nominal trajectory
%
% Inputs:
%   x_nom -> nominal state
%   u_nom -> nominal control
%
% Outputs:
%   A -> linearized A matrix
%   B -> linearized B matrix
%   C -> linearized C matrix
%   D -> linearized D matrix
%
% Author: Thomas Dunnington
% Modified: 12/2/2024

% Constant parameters
L = 0.5;

% State variables
xi_g = x_nom(1); eta_g = x_nom(2); theta_g = x_nom(3); xi_a = x_nom(4); eta_a = x_nom(5); theta_a = x_nom(6);

% Control inputs
v_g = u_nom(1); phi_g = u_nom(2); v_a = u_nom(3); w_a = u_nom(4);

% Jacobians
A = [0, 0, -v_g*sin(theta_g), 0, 0, 0;
    0, 0, v_g*cos(theta_g), 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, -v_a*sin(theta_a);
    0, 0, 0, 0, 0, v_a*cos(theta_a);
    0, 0, 0, 0, 0, 0];

B = [cos(theta_g), 0, 0, 0;
    sin(theta_g), 0, 0, 0;
    1/L * tan(phi_g), v_g/L * (sec(phi_g))^2, 0, 0;
    0, 0, cos(theta_a), 0;
    0, 0, sin(theta_a), 0;
    0, 0, 0, 1];


% ADD PROPER JACOBIANS HERE, THESE ARE PLACEHOLDERS
C = [(eta_a-eta_g)/((eta_a-eta_g)^2+(xi_a-xi_g)^2), -(xi_a-xi_g)/((eta_a-eta_g)^2+(xi_a-xi_g)^2), -1, -(eta_a-eta_g)/((eta_a-eta_g)^2+(xi_a-xi_g)^2), (xi_a-xi_g)/((eta_a-eta_g)^2+(xi_a-xi_g)^2), 0;
    (xi_g-xi_a)/(sqrt((xi_g-xi_a)^2+(eta_g-eta_a)^2)), (eta_g-eta_a)/(sqrt((xi_g-xi_a)^2+(eta_g-eta_a)^2)),0,-(xi_g-xi_a)/(sqrt((xi_g-xi_a)^2+(eta_g-eta_a)^2)), -(eta_g-eta_a)/(sqrt((xi_g-xi_a)^2+(eta_g-eta_a)^2)), 0;
    -(eta_g-eta_a)/((eta_g-eta_a)^2+(xi_g-xi_a)^2), (xi_g-xi_a)/((eta_g-eta_a)^2+(xi_g-xi_a)^2), 0, (eta_g-eta_a)/((eta_g-eta_a)^2+(xi_g-xi_a)^2), -(xi_g-xi_a)/((eta_g-eta_a)^2+(xi_g-xi_a)^2), -1;
    0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 1, 0];

D = [0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0];


end


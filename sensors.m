function y = sensors(x)
%SENSORS Given a state vector x, this function calculates a measurement
%with no measurement noise
%
% Inputs:
%   x -> full state vector (xi_g, eta_g, theta_g, xi_a, eta_a, theta_a)
%
% Outputs: 
%   y -> measurement vector (azimuth relative to UGV, relative distance,
%   azimuth realtive to UAV, UAV easting, UAV northing)
%
% Author: Thomas Dunnington
% Modified: 12/2/2024

% Azimuth measurements
xi_g = x(1); eta_g =x(2); theta_g = x(3); xi_a = x(4); eta_a = x(5); theta_a = x(6);
azimuth_ugv = atan((eta_a - eta_g) / (xi_a - xi_g)) - theta_g;
azimuth_uav = atan((eta_g - eta_a) / (xi_g - xi_a)) - theta_a;

% Range measurement
range = sqrt((xi_g - xi_a)^2 + (eta_g - eta_a)^2);

% Full measurement vector
y = [azimuth_ugv; range; azimuth_uav; xi_a; eta_a];


end


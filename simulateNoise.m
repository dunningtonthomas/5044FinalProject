function [time,x_noise_mat,y_noise_mat] = simulateNoise(x0,u0,Q,R,dt,n_ind)
rng(100);
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
Sv_process = chol(Q,'lower');
x_noise_mat = x0';
time = 0;
for i=1:n_ind
    % Define time step
    ZOH_time = dt*[i-1 i];
    % Generate AWGN for the time step
    q = randn(6,1);
    w =  Sv_process*q;
    % Set Initial Conditions
    x_init = x_noise_mat(end,:);
    % Function Handle for the nominal trajectory using the full nonlinear equations with noise
    eomFunc = @(t, x)coopEOM(t, x, u0, w);
    % Simulate over time step
    [TOUT, x_noise_step] = ode45(eomFunc, ZOH_time, x_init, options);
    x_noise_mat = [x_noise_mat;x_noise_step(end,:)];
    time = [time;TOUT(end)];
end

% Angle wrapping
x_noise_mat(:,3) = mod(x_noise_mat(:,3) + pi, 2*pi) - pi;
x_noise_mat(:,6) = mod(x_noise_mat(:,6) + pi, 2*pi) - pi;

% Calculate the measurements from the sensor model
y_noise_mat = zeros(length(time), 5);
for i = 1:length(time)
    y_noise_mat(i,:) = sensors(x_noise_mat(i,:))';
end
% Add Measurement Noise 
Sv_measurement = chol(R,'lower');
q = randn(5,length(time));
y_noise_mat =  y_noise_mat + (Sv_measurement*q)';
end


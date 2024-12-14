function [time,x_noise_mat,y_noise_mat] = simulateNoise(x0,u0,Q,R,dt,n_ind)
    % LKF This function is a Linearized Kalman Filter that returns the estimated state, measurements and covariance
    %
    % Inputs: 
    %   x0          -> Initial State without noise
    %   u0          -> Controls
    %   Q           -> Process Noise Covariance
    %   R           -> Measurement Noise Covariance
    %   dt          -> time step
    %   n_ind       -> Number of Data Points after the Initial State (if you input 1000 it will return 1000 new points after the initial state so the final output size is 1001)
    % Outputs:
    %   time        -> Simulation Time Vector 
    %   x_noise     -> Simulated State over Time with Process Noise
    %   y_noise     -> Estimated Standard Deviation at all time steps
    %
    % Author: Owen Craig
    % Modified: 12/3/2024
    % Define Constant RNG Seed and ODE45 Tol
    %rng(100);
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    
    % Simulate Process Noise
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
    
    % Calculate the measurements from the sensor model
    y_noise_mat = zeros(length(time)-1, 5);
    for i = 2:length(time)
        y_noise_mat(i-1,:) = sensors(x_noise_mat(i,:))';
    end
    % Add Measurement Noise 
    Sv_measurement = chol(R,'lower');
    q = randn(5,length(time)-1);
    y_noise_mat =  y_noise_mat + (Sv_measurement*q)';
end


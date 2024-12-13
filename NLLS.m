function [xhat_final,P_final] = NLLS(u, y_actual, num_measurements, Q, R, dt)
%NLLS Use nonlinear least squares for a warm start for the EKF to estimate
%the initial state and covariance
%
% Inputs: 
%   u -> constant control input
%   y_actual -> simulated measurement data, p x num_meas matrix
%   num_measurements -> amount of measurements to use from y_actual
%   Q -> process noise covariance matrix
%   R -> sensor noise covariance matrix
%   dt -> time step size
% Outputs:
%   xhat_meas -> next full state estimate
%   P_meas -> next covariance
%
% Author: Thomas Dunnington
% Modified: 12/12/2024

    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
    % Nonlinear least squares warm start
    
    % Create sensor noise matrix
    R_cells = repmat({R}, 1, num_measurements);
    Rbig = blkdiag(R_cells{:});
    
    % Get 10 measurements of data and stack the data
    y_meas = [];
    for i = 1:num_measurements
        y_meas = [y_meas; y_actual(i,:)'];
    end
    
    % Initial guess
    xhat_old = [0; 0; 0; y_meas(4); y_meas(5); 0]; % Rough initial guess, maybe use the first measurement
    
    % Gauss Newton optimization
    max_iterations = 200;
    alpha_splits = 10;
    iteration = 0;
    time_vector = 0:dt:(num_measurements*dt);
    
    while iteration < max_iterations
        % Get predicted measurement vector
        y_pred = [];
        Hbig = [];
        xinit = xhat_old;
        [~, xpred] = ode45(@(t, x)coopEOM(t, x, u, zeros(length(xinit),1)), time_vector, xinit, options);
    
        for i = 1:num_measurements
            % Get predicted measurement
            y_pred = [y_pred; sensors(xpred(i+1, :)')];
    
            % Measurement jacobian
            [~, ~, H, ~] = linearize(xpred(i+1, :)', zeros(4,1));
            Hbig = [Hbig; H];
        end
    
        % Calculate the cost function
        Jcurr = (y_meas - y_pred)' * inv(Rbig) * (y_meas - y_pred);
        dxhat = inv(Hbig' * inv(Rbig) * Hbig) * Hbig' * inv(Rbig) * (y_meas - y_pred);
    
        % Make alpha smaller so we don't diverge
        alpha = 1;
        alpha_iter = 1;
        while alpha_iter < alpha_splits
            % Proposed update to the estimate
            xhat_new = xhat_old + alpha .* dxhat;
            
            % Get new predicted measurements
            yhat_new = [];
            xinit = xhat_new;
            [~, xpred] = ode45(@(t, x)coopEOM(t, x, u, zeros(length(xinit),1)), time_vector, xinit, options);
    
            for i = 1:num_measurements
                % Predicted sensor measurements
                yhat_new = [yhat_new; sensors(xpred(i+1, :)')];
            end
    
            % Calculate the cost function
            Jnew = (y_meas - yhat_new)' * inv(Rbig) * (y_meas - yhat_new);
    
            % Take best cost
            if Jnew > Jcurr
                alpha = alpha / 2; % Shrink size
            else
                break;
            end
    
        end
    
        % Output the iteration and the cost function
        fprintf('Iteration: %d, Cost: %.6f\n', iteration, Jcurr);
    
        % Update the best guess
        xhat_old = xhat_new;
        iteration = iteration + 1;
    end
    
    % Final estimate
    xhat_final = xhat_old;
    
    % Calculate the final covariance
    [~, ~, H, ~] = linearize(xhat_final, zeros(4,1));
    P_final = inv(H'*H);
end


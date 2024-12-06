function [xhat_final,P_final] = NLLS(u,y_actual,Q,R,dt)
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    % Nonlinear least squares warm start
    num_measurements = 10;
    
    % Create sensor noise matrix
    R_cells = repmat({R}, 1, num_measurements);
    Rbig = blkdiag(R_cells{:});
    
    % Get 10 measurements of data and stack the data
    y_meas = [];
    for i = 1:num_measurements
        y_meas = [y_meas; y_actual(:,i+1)];
    end
    
    % Initial guess
    xhat_old = [0; 0; 0; 0; 0; 0]; % Rough initial guess, maybe use the first measurement
    
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


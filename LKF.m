function [x_est,sigma]= LKF(x_nom,u_nom,y_nom,y_actual,u_actual,Q,R,dt)
    % LKF This function is a Linearized Kalman Filter that returns the estimated state, measurements and covariance
    %
    % Inputs: 
    %   x_nom       -> Nominal State Trajectory
    %   u_nom       -> Nominal Control Inputs
    %   y_nom       -> Nominal Measured States
    %   y_actual    -> Actual measured state at all times
    %   u_actual    -> Actual control inputs at all times
    %   Q           -> Process Noise Covariance
    %   R           -> Measurement Noise Covariance
    %   dt          -> time step
    % Outputs:
    %   x_est       -> Estimated state for all time steps
    %   y_update    -> Estimates measurement for all time steps
    %   sigma       -> Estimated Standard Deviation at all time steps
    %
    % Author: Owen Craig
    % Modified: 12/3/2024
    n = length(y_actual);
    H_mat = [];
    Gamma = eye(size(x_nom,1),size(x_nom,1));
    for i = 1:n
        [A, B, C, D] = linearize(x_nom(:,i), u_nom(:,i));
        [~,G,~,H] = eulerDiscretize(A,B,C,D,Gamma,dt);
        H_mat = [H_mat;H];
    end
    % Initialize the filter using a batch lls to warm start the filter @k = 0
    % R_mat = kron(eye(n), R);
    % y_reshaped = reshape(y_actual', [], 1);
    % xls  = (H_mat'*R_mat^-1*H_mat)^-1*H_mat'*R_mat^-1*y_reshaped;
    % Pls = (H_mat'*R_mat^-1*H_mat)^-1;
    % P_update = Pls;
    % dx_update = xls-x_nom(:,1);
    % dx_update(3) = mod(dx_update(3) + pi, 2*pi) - pi;
    % dx_update(6) = mod(dx_update(6) + pi, 2*pi) - pi;
    % du = zeros(size(G,2),1);
    % dy_update = y_actual-y_nom;
    % x_est(:,1) = dx_update+x_nom(:,1);
    % sigma(:,1) = sqrt(diag(P_update));
     x_est(:,1) = x_nom(:,1);
    
     P_update = diag(ones(6,1)*100000);
     sigma(:,1) = sqrt(diag(P_update));
     dx_update = x_est(:,1)-x_nom(:,1);
     dx_update(3) = mod(dx_update(3) + pi, 2*pi) - pi;
     dx_update(6) = mod(dx_update(6) + pi, 2*pi) - pi;
     du = zeros(size(G,2),1);
     dy_update = y_actual-y_nom;
     
    % Use the lineraized Kalman Filter to estimate the state
    for i =1:n
        % Pure Predition uses F_k and G_k so calculate @ k
        [A, B, C, D] = linearize(x_nom(:,i), u_nom(:,i));
        [F,G,~,~] = eulerDiscretize(A,B,C,D,Gamma,dt);
        % Measurement Update uses H_K+1 so calculate H @ k+1
        [A, B, C, D] = linearize(x_nom(:,i+1), u_nom(:,i+1));
        [~,~,Omega,H] = eulerDiscretize(A,B,C,D,Gamma,dt);
        % Prediction Step
        % Estimate state perturbation
        dx = F*dx_update+G*du;
        P = F*P_update*F'+Omega*Q*Omega';
        du = u_actual(:,i+1) - u_nom(:,i+1);
        % Angle Wrapping
        dx(3,:) = mod(dx(3,:) + pi, 2*pi) - pi;
        dx(6,:) = mod(dx(6,:) + pi, 2*pi) - pi;
        % Measurement Update Step
        K = P*H'*(H*P*H'+R)^-1;
        dx_update = dx +K*(dy_update(:,i)-H*dx);
        P_update = (eye(size(P))-K*H)*P;
        % Angle Wrapping
        dx_update(3,:) = mod(dx_update(3,:) + pi, 2*pi) - pi;
        dx_update(6,:) = mod(dx_update(6,:) + pi, 2*pi) - pi;

        % Save Estimates
        x_est(:,i+1) = dx_update+x_nom(:,i);
        sigma(:,i+1) = sqrt(diag(P_update));
    end
    x_est(3,:) = mod(x_est(3,:) + pi, 2*pi) - pi;
    x_est(6,:) = mod(x_est(6,:) + pi, 2*pi) - pi;
end
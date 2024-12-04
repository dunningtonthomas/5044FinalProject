function [x_est,sigma]= LKF(x_nom,u_nom,y_nom,y_actual,u_actual,Q,R,dt)
    % LKF This function is a Linearized Kalman Filter that returns the estimated state, measurements and covariance
    %
    % Inputs: 
    %   x_nom       -> DT State Space
    %   u_nom       -> Sensor Noise Covariance
    %   y_nom       -> Process Noise Covariance
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
    n = length(x_nom);
    H_mat = [];
    Gamma = eye(size(x_nom,1),size(x_nom,1));
    for i = 1:n
        [A, B, C, D] = linearize(x_nom(:,i), u_nom(:,i));
        [F,G,Omega,H] = eulerDiscretize(A,B,C,D,Gamma,dt);
        H_mat = [H_mat;H];
    end
    % Initialize the filter using a batch lls to warm start the filter
    R_mat = kron(eye(n), R);
    y_reshaped = reshape(y_actual', [], 1);
    xls  = (H_mat'*R_mat^-1*H_mat)^-1*H_mat'*R_mat^-1*y_reshaped;
    Pls = (H_mat'*R_mat^-1*H_mat)^-1;
    P_update = Pls;
    dx_update = xls-x_nom(:,1);
    dx_update(3) = mod(dx_update(3) + pi, 2*pi) - pi;
    dx_update(6) = mod(dx_update(6) + pi, 2*pi) - pi;
    du = zeros(size(G,2),1);
    dy_update = y_actual-y_nom;
    % Use the lineraized Kalman Filter to estimate the state
    for i =1:n
        [A, B, C, D] = linearize(x_nom(:,i), u_nom(:,i));
        [F,G,Omega,H] = eulerDiscretize(A,B,C,D,Gamma,dt);
        % Prediction Step
        % Estimate state perturbation
        dx = F*dx_update+G*du;
        P = F*P_update*F'+Omega*Q*Omega;
        du = u_actual(:,i) - u_nom(:,i);
        % Measurement Update Step
        K = P*H'*(H*P*H'+R)^-1;
        dx_update = dx +K*(dy_update(:,i)-H*dx);
        P_update = (eye(size(P))-K*H)*P;
        x_est(:,i) = dx_update+x_nom(:,i);
        sigma(:,i) = sqrt(diag(P_update));
    end

end
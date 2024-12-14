function [x_est,sigma,innovation_vec,S_vec]= LKF(x_nom,u_nom,y_nom,y_actual,u_actual,Q,R,dt)
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
    Gamma = eye(size(x_nom,1),size(x_nom,1));
    % Initialize Filter
    x_est(:,1) = x_nom(:,1);
    P_update = Q*10;
    % Initialize Output Variables
    S_vec = {};
    innovation_vec = [];
    sigma(:,1) = sqrt(diag(P_update));
    dx_update = x_est(:,1)-x_nom(:,1);
    du = zeros(size(u_nom,1),1);
    dy_update = y_actual-y_nom;
    dy_update(1) = mod(dy_update(1) + pi, 2*pi) - pi;
    dy_update(3) = mod(dy_update(3) + pi, 2*pi) - pi;
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
        dx = F*dx_update;
        P = F*P_update*F'+Omega*Q*Omega';
        % Measurement Update Step
        S = H*P*H'+R;
        K = P*H'*(S)^-1;
        innovation = dy_update(:,i)-H*dx;
        % Angle wrap the innovation
        innovation(1) = mod(innovation(1) + pi, 2*pi) - pi;
        innovation(3) = mod(innovation(3) + pi, 2*pi) - pi;
        dx_update = dx +K*(innovation);
        P_update = (eye(size(P))-K*H)*P;

        x_est(:,i+1) = dx_update+x_nom(:,i+1);
        S_vec{i} = S;
        innovation_vec(:,i) = innovation;
        sigma(:,i+1) = sqrt(diag(P_update));
    end
end
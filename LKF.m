function [x_est,sigma,innovation_vec,S_vec,P_vec]= LKF(x_nom,u_nom,y_nom,y_actual,u_actual,Q,R,dt)
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
    %% Define constants
    Gamma = eye(size(x_nom,1),size(x_nom,1));
    du = zeros(size(u_nom,1),1);
    dy_update = y_actual-y_nom;
    dy_update(1) = mod(dy_update(1) + pi, 2*pi) - pi;
    dy_update(3) = mod(dy_update(3) + pi, 2*pi) - pi;
    n = length(y_actual);
    % Initialize Output Variables
    S_vec = {};
    innovation_vec = [];
    
    
    % Initialize Filter
    % x_est(:,1) = x_nom(:,1);
    % P_update = diag([10, 10, pi/10, 5, 5, pi/10]);
    
    % WARM START THE FILTER USING BATCH LLS FOR THE FIRST 4 MEASUREMENTS
    H_mat = [];
    samples = 6;
    for i = 1:samples
        [A, B, C, D] = linearize(x_nom(:,i), u_nom(:,i));
        [~,G,~,H] = eulerDiscretize(A,B,C,D,Gamma,dt);
        H_mat = [H_mat;H];
    end
    R_mat = kron(eye(samples), R);
    y_reshaped = dy_update(:,1:samples);
    y_reshaped = y_reshaped(:);
    dxls  = (H_mat'*R_mat^-1*H_mat)^-1*H_mat'*R_mat^-1*y_reshaped;
    Pls = (H_mat'*R_mat^-1*H_mat)^-1;
    P_update = Pls;
    P_update = diag([10, 10, pi/10, 1, 1, pi/10]);
    %dx_update = dxls;
    dx_update = zeros(size(dxls));
    % dx_update(3) = mod(dx_update(3) + pi, 2*pi) - pi;
    % dx_update(6) = mod(dx_update(6) + pi, 2*pi) - pi;
    x_est(:,1) = dx_update+x_nom(:,1);
    sigma(:,1) = sqrt(diag(P_update));
    P_vec{1} = P_update;
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
        P_vec{i} = P_update;
    end
end
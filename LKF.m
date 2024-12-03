function [x_update,P_update]= LKF(x_nom,u_nom,y_nom,Q,R,dt)
    % LKF This function is a Linearized Kalman Filter that returns the estimated state, measurements and covariance
    %
    % Inputs: 
    %   F,G,H,M     -> DT State Space
    %   R           -> Sensor Noise Covariance
    %   Q           -> Process Noise Covariance
    %   x_meas      -> Last state estimate after the measurement update at time k
    %   y           -> Measured state at time k+1
    %   P           -> Last covariance update at time k
    %   u           -> Control input at time k
    % Outputs:
    %   x_update   -> Estimated state perturbation from the nominal state at time k+1
    %   y_update   -> Estimates measurement perturbation from nominal measurement at time k+1
    %   P_update    -> Estimated state covariance at time k+1
    %
    % Author: Owen Craig
    % Modified: 12/3/2024
    [A, B, C, D] = linearize(x_nom, u_nom);
    Gamma = eye(size(A));
    [F,G,Omega,H] = eulerDiscretize(A,B,C,D,Gamma,dt);
    % Prediction Step
    % Estimate state perturbation
    dx = F*dx_update+G*du;
    P = F*P_update*F'+Omega*Q*Omega;
    du = u_actual - u_nom;
    % Measurement Update Step
    dx_update = dx +K*(dy_update-H*dx);
    P_update = (eye(P)-K*H)*P;
    K = P*H'*(H*P*H'+R)^-1;
    dy_update = y_actual-y_nom;
end
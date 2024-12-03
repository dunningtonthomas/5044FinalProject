function [F,G,Omega,H] = eulerDiscretize(A,B,C,D,Gamma,dt)
    %eilerDiscretize This function uses eulers method to approximate the DT State Space.
    %
    % Inputs: 
    %   A,B,C,D     -> CT Linearized State Space
    %   Gamma       -> Process noise map to linearized CT state space
    %   dt          -> Time step
    % Outputs:
    %   F,G,H       -> DT Linearized State Space
    %   Omega       -> Process noise map to linearized DT state space
    % Author: Owen Craig
    % Modified: 12/3/2024

    % Using Eulers method convert CT perturbation model into DT model (Eq. from lecture 31)
    F = eye(size(A))+dt*A;
    G = dt*B;
    Omega = dt*Gamma;
    H = C;
end


function [F, G, H, M] = discretize(A, B, C, D, dt)
%LINEARIZE This function finds the DT linearized model given a CT linear
%model
%
% Inputs:
%   A, B, C, and D matrices of the CT linear model
%   dt -> ZOH approximation time step
%
% Outputs:
%   F, G, H, and M matrices of the DT linear model
%
% Author: Thomas Dunnington
% Modified: 12/2/2024

% F and G matrices
Ahat = [A, B;
    zeros(length(B(1,:)), length(A(:,1)) + length(B(1,:)))];
expm_hat = expm(Ahat.*dt);
F = expm(A.*dt);
G = expm_hat(1:length(A(:,1)), length(A(:,1))+1:end);

% H and M
H = C;
M = D;

end


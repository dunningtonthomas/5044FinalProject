function sigma_points = sigmaPoints(xhat, P, lambda2)
% Generate sigma points for the given state and covariance
S = chol(P, 'lower');
sigma_points = [xhat, xhat + lambda2 * S, xhat - lambda2 * S];
end
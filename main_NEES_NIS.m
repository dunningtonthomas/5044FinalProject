%% Clean
close all; clear; clc;

%% LFK
% Setup EOM
% Ode45 Constants
Data = load('cooplocalization_finalproj_KFdata.mat');
dt = 0.1;
tspan = [0 100];
t_nom = (dt:dt:tspan(2))';
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

% Nominal values
x_ugv = [10; 0; pi/2];
x_uav = [-60; 0; -pi/2];
u_ugv = [2; -pi/18];
u_uav = [12; pi/25];

x_nom = [x_ugv; x_uav];
u_nom = [u_ugv; u_uav];

% Simulate Nominal Nonlinear Trajectoy without Noise
w = zeros(6,1);
eomFunc = @(t, x)coopEOM(t, x, u_nom, w);
x_init = x_nom;
t_nom = (0:dt:tspan(2))';
[~, x_nom_mat] = ode45(eomFunc, t_nom, x_init, options);
u_nom_mat = ones(length(t_nom), 4) .* u_nom';

% Calculate the measurements from the sensor model
y_nom_mat = zeros(length(t_nom)-1, 5);
for i = 2:length(t_nom)
    y_nom_mat(i-1,:) = sensors(x_nom_mat(i,:))';
end

% Simulate Nonlinear Trajectory with Process Noise
Q_true = Data.Qtrue;
R_true = Data.Rtrue;

% Apply Linearized Kalman Filter
Q_tune = Q_true;
Q_tune(1,1) = Q_tune(1,1)*.1;
Q_tune(2,2) = Q_tune(2,2)*.1;
Q_tune(3,3) = Q_tune(3,3)*1000000;
Q_tune(4,4) = Q_tune(4,4)*0.1;
Q_tune(5,5) = Q_tune(5,5)*0.1;
Q_tune(6,6) = Q_tune(6,6)*0.00001;


Q_tune(1,2) = -6/1000000;
Q_tune(2,1) = Q_tune(1,2);

Q_tune(1,3) = 1/10000;
Q_tune(3,1) = Q_tune(1,3);

Q_tune(2,3) = 1/100000;
Q_tune(3,2) = Q_tune(2,3);

Q_tune(1,4) = 1/100000;
Q_tune(4,1) = Q_tune(1,4);

Q_tune(1,5) = 1/100000;
Q_tune(5,1) = Q_tune(1,5);

Q_tune(2,4) = 1/100000;
Q_tune(4,2) = Q_tune(2,4);

Q_tune(2,5) = 1/100000;
Q_tune(5,2) = Q_tune(2,5);

Q_tune = Q_tune.*200000;



% Monte Carlo Simulations for NEES and NIS
Nsim = 50; % Number of Monte Carlo simulations
Nstate = size(x_nom, 1);
Nmeas = size(y_nom_mat, 2);

% Preallocate for NEES and NIS results
nees_values = zeros(Nsim, length(t_nom)-1);
nis_values = zeros(Nsim, length(t_nom)-1);

% Initialize P 
 P_init = diag([1, 1, pi/10, 1, 1, pi/10]);


for sim_idx = 1:Nsim
    Sv = chol(P_init,'lower');
    q = randn(size(P_init,1),1);
    x_sim_init = x_init+Sv*q;

    % Simulate noisy trajectory
    [~, x_noisy, y_noisy] = simulateNoise(x_sim_init, u_nom, Q_true, R_true, dt, 1000);

    % Apply Linearized Kalman Filter
    [x_LKF, sigma, innovation_vec, S_vec,P_vec] = LKF(x_nom_mat', u_nom_mat', y_nom_mat', y_noisy', u_nom_mat', Q_tune, R_true, dt);

    % Compute NEES and NIS for this run
    for k = 2:length(t_nom)-1
        % State estimation error
        e_k =  x_LKF(:, k)-x_noisy(k, :)';
        e_k(3) = mod(e_k(3)+pi,2*pi)-pi;
        e_k(6) = mod(e_k(6)+pi,2*pi)-pi;
        P_k = P_vec{k};

        % NEES (normalized state error)
        nees_values(sim_idx, k) = e_k' * (P_k \ e_k); % e_k' * inv(P_k) * e_k

        % NIS (normalized innovation error)
        innovation = innovation_vec(:,k);
        innovation(1) = mod(innovation(1)+pi,2*pi)-pi;
        innovation(3) = mod(innovation(3)+pi,2*pi)-pi;
        S = S_vec{k};
        nis_values(sim_idx, k) = innovation' * (S \ innovation); % innov' * inv(S_k) * innov
    end
end

% Chi-Square Test Bounds
alpha = 0.05; % Significance level
r1_NEES = chi2inv(alpha/2, Nsim*Nstate) / Nsim;
r2_NEES = chi2inv(1-alpha/2, Nsim*Nstate) / Nsim;

r1_NIS = chi2inv(alpha/2, Nsim*Nmeas) / Nsim;
r2_NIS = chi2inv(1-alpha/2, Nsim*Nmeas) / Nsim;

% Average NEES and NIS across all runs
mean_nees = mean(nees_values, 1); % Time-averaged NEES
mean_nis = mean(nis_values, 1);   % Time-averaged NIS


% Plotting
figure;
plot(t_nom(2:end), mean_nees, 'ro', 'LineWidth', 1.5);
hold on;
yline(r1_NEES, 'r--', 'LineWidth', 1.2);
yline(r2_NEES, 'r--', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('NEES');
legend('Mean NEES', '\chi^2 Lower Bound', '\chi^2 Upper Bound');
title('NEES Chi-Square Test for LKF');
grid on;

figure
plot(t_nom(2:end), mean_nis, 'bo', 'LineWidth', 1.5);
hold on;
yline(r1_NIS, 'b--', 'LineWidth', 1.2);
yline(r2_NIS, 'b--', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('NIS');
legend('Mean NIS', '\chi^2 Lower Bound', '\chi^2 Upper Bound');
title('NIS Chi-Square Test for LKF');
grid on;

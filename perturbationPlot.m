function perturbationPlot(TOUT, XOUT, XNOM, linespec)
%PERTURBATIONPLOT Plot the perturbations of each state given a nominal
%trajectory and the actual trajectory
%states over time
% 
% Inputs:
%   TOUT -> time vector
%   XOUT -> length(TOUT) by 6 state matrix, each column is a state
%   XNOM -> nominal trajectory matrix
%   linespec -> linestyle of the plot
%
% Output:
%   plot 1 -> 6x1 subplot of the state perturbations
%
% Author: Thomas Dunnington
% Modified: 12/3/2024

% Calculate the perturbations
XPERT = XOUT - XNOM;


figure();
subplot(6,1,1)
plot(TOUT, XPERT, )



end


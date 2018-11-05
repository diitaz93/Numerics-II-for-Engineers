% Assignment 01, Exercise 4c) Timing as n tends to infinity
% by Shane Alpert, Juan Sebastian Diaz, Bo Yan
% DESCRIPTION: This script plots the time that it takes for our method to
% compute the solution for the linear system Kn*x=b for several different
% sizes.
clc;
close all
n = [ logspace(1,7,7)]; % Vector with logarithmic scaling of size
time = zeros(1,7); % Initialization of time vector

% 'Warming up' the 'just in time' compiler
b = randn(n(1),1);
a01e04invKn(b);

% Loop running over sizes
for i=1:7
    b = randn(n(i),1); % b initialized randomly
    tic;
    a01e04invKn(b); % Solves for x in Kn*x=b
    time(i) = toc; % Save timing in time vector
end
% Logarthimic plotting
loglog(n,time,'+-');

title('Execution time vs. size of matrix')
xlabel('Size (n)')
ylabel('Time (s)')
%% Discussion
% When n approaches infinity, the time increases linearly with n. For small
% n, times are not linear, probably due to 'just in time' or compiler
% issues.
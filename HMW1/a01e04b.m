% Assignment 01, Exercise 4b) Comparision of times and residuals
% by Shane Alpert, Juan Sebastian Diaz, Bo Yan
% DESCRIPTION: This script compares the execution times and residuals of 
% our method to solve the linear system Kn*n=b and three matlab functions.
% The times and residuals are plotted in two figures.
clc;
close all
n = [ 10 100 1000]; % Vector of different sizes of matrices Kn
l = length(n); % Represents the number of sizes used in for loop
% Initialization of Norm matrix. Rows represent the different methods and
% columns the respective size
N = zeros(4,l); 
% Initialization of 4 time vectors, one for each method
time_our = zeros(1,l);
time_1 = zeros(1,l);
time_2 = zeros(1,l);
time_3 = zeros(1,l);

% Loop running over the sizes
for i=1:l
    b = randn(n(i),1); % b initialized as random vector of n(i) entries
    %% Our function
    tic;
    % (Kn is computed inside the function)
    x = a01e04invKn(b); % x computation
    time_our(i) = toc;
    N(1,i) = norm(full(a01e03sparse(n(i)))*x-b); % Norm computation
    %% Function 1: Kn\b
    tic;
    Kn = a01e03sparse(n(i)); % Kn computation
    x = Kn\b; % x computation
    time_1(i) = toc;
    N(2,i) = norm(full(Kn)*x-b); % Norm computation
    %% Function 2: full(Kn)\b
    tic;
    Kn = a01e03sparse(n(i)); % Kn computation
    x = full(Kn)\b; % x computation
    time_2(i) = toc;
    N(3,i) = norm(full(Kn)*x-b); % Norm computation
    %% Function 3: inv(full(Kn))*b
    tic;
    Kn = a01e03sparse(n(i)); % Kn computation
    x = inv(full(Kn))*b; % x computation
    time_3(i) = toc;
    N(4,i) = norm(full(Kn)*x-b); % Norm computation
end
%% Time vs. dimension logarithmic plotting
loglog(n,time_our,'b')
title('Execution time for different methods')
xlabel('Dimension of vector (n)')
ylabel('Time (s)')
hold on
loglog(n,time_1,'m')
loglog(n,time_2,'r')
loglog(n,time_3,'g')
legend('Our function','Kn\b','full(Kn)\b','inv(full(Kn))*b','location'...
    ,'northwest')
%% Norm vs. Dimension logarithmic plotting
figure;
loglog(n,N(1,:),'b')
title('Norm 2 of residual for different methods')
xlabel('Dimension of vector (n)')
ylabel('Norm 2 of residual')
hold on
loglog(n,N(2,:),'m')
loglog(n,N(4,:),'r')
loglog(n,N(3,:),'g')
legend('Our function','Kn\b','full(Kn)\b','inv(full(Kn))*b','location'...
    ,'northwest')
%% Discussion
% As it can be seen in the graphics (run several times) our method beats
% the MATLAB function for all n and has the lowest norm of the residual as
% well. This is because our function is optimized for the specific matrix
% to be used (Kn), but MATLAB functions are designed to solve the problem 
% for any type of matrix.
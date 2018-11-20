function [ error ] = a03ex03solve( )
% Assignment 03, Exercise 03c, by Shane Alpert (404579), 
% Juan Sebastian Diaz (405385), Bo Yan (403787)
% 
% Solves the discretized BVP problem for different values of p, which 
% determines the number of grid points and the mesh size. For each
% different p, the error between the approximation(numerical solution) and the restricted exact
% solution in the maximum norm is determined. 

p = 1:15;   % Creates array of integer p values from 1 to 15
N = 2.^p-1; % number of internal grid points
h = 1./(N+1);% mesh size
error = zeros(1,15);% creates array for errors and initialize it

for i = 1:length(p) % solves discretized BVP for different mesh sizes
    [ xh, Lh, fh ] = a03ex03getBVP(p(i)); % discretize the BVP problem
    uh = Lh \ fh; % Obtain numerical solution
    u = 1 + 4*xh.^2 - 3*xh.^3; %Obtain exact values from exact solution
    % u is a vector with the same dimension as uh has
    error(i) = max(abs((uh - u)));% store the error
end

loglog(h,error);
xlabel('h');
ylabel('error');
K=polyfit(log(h),log(error),1)% Fit the polynomial with the data.
%K[1] is the order of convergence.
%10^{K[2]} is the coefficient
%i.e.error(p) is equivalent to 10^{K[2]}*h^K[1].
%Run the program we can get the answer. K[1]=2.
%Thus, the order of convergence is 2.
end


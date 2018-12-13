function [Lh, fh] = a07e04getPDE(p,beta)
%% Assignment 06, Exercise 02, by Shane Alpert, Juan Sebastian Diaz, Bo Yan
%a07e03getPDE Calculates the Lh matrix and fh vector of the two-point
%boundary problem with Robin boundary conditions. We calculate Lh to be of
%the form Lh = A*S+B*T+C
N = 2^p-1;
h = 2/(N+1);
%% A(x)
i = 1:N;
a = 2-(h*i-1).^2;
r = 1:N;
c = 1:N;
A = sparse(r,c,a,N,N);

%% S
r = [1:N 1:N-1 2:N];
c = [1:N 2:N 1:N-1];
s = [[-2*ones(1,N-1) -1] ones(1,2*N-2)];
S = sparse(r,c,s,N,N);

%% B(x)
b = -h*i+1;
r = 1:N;
c = 1:N;
B = sparse(r,c,b,N,N);

%% T
r = [1:N 1:N-1 2:N];
c = [1:N 2:N 1:N-1];
t = [[zeros(1,N-1) 1] ones(1,N-1) -1*ones(1,N-1)];
T = sparse(r,c,t,N,N);

%% C
C = 16*speye(N);

%% Lh
Lh = 1/(h^2)*A*S+1/(2*h)*B*T+C;

%% fh

fh = zeros(N,1);
fh(1) = (-h^2+3*h+2)/4*(h^2);
fh(end) = (2*(N*h)^2-4*N*h-2+N*h^2-h)*beta/(2*h);
end

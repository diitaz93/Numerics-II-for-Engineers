function [Lh,fh] =  a06e02getPDE(N,f,g)
% a06e02getPDE solves a PDE numerically given the size of the grid N, and
% the function handles for the inhomogeinity function f and the initial
% condition g. It returns the matrix Lh and the function f discretized in
% the vector fh.
% Assignment 06, Exercise 03, by Shane Alpert (404579), 
% Juan Sebastian Diaz (405385), Bo Yan (403787) 

%% Grid parameters
h=1/(N+1);
%% Lh
% Create T matrix
r = [1:N 1:N-1 2:N];
c = [1:N 2:N 1:N-1];
v = [[-3 -4*ones(1,N-2) -3] 1*ones(1,2*N-2)];
T = sparse(r,c,v,N,N);
%
I = eye(N);
% Allocate memory for a sparse matrix of size N^2xN^2
Lh = spalloc(N^2, N^2,(N-2)*nnz(T)+2*N*(N-1)+2*nnz(T+I));
% Build matrix
U = [T+I I];
M = [I T I];
B = [I T+I];
% Set the upper ando lower rows manually
Lh(1:N,1:2*N) = U;
Lh((N^2-N)+1:N^2 , (N^2-2*N)+1:N^2) = B;
% Fill the middle of the matrix
for j=1:N-2
    Lh(j*N+1:(j+1)*N, (j-1)*N+1:(j+2)*N) = M;
end
Lh=-Lh/(h^2);
 
%% fh=Rhf+G
% Rhf (funtion evaluated in the grid
[X,Y] = meshgrid(h:h:1-h);
Rhf= f(X,Y);
% G
% Column 1, all rows, i=1
x1 = zeros(N,1);
y1 = (1:N)'*h;

g1 = g(x1,y1);
G1 = [g1 zeros(N,N-1)];

% Row 1, all columns,j=1
x2 = [1:N]*h;
y2 = zeros(1,N);

g2 = g(x2,y2);
G2 = [g2; zeros(N-1,N)];

% Column N, all rows, i=N
x3 = ones(N,1);
y3 = [1:N]'*h;

g3 = g(x3,y3);
G3 = [zeros(N,N-1) g3];

% Row N, all columns, j=N
x4 = [1:N]*h;
y4 = ones(1,N);

g4 = g(x4,y4);
G4 = [zeros(N-1,N); g4];
%G
G=(G1 + G2 + G3 + G4)/h;
%fh
fh = Rhf+ G;
fh = fh';
fh = fh(:);

end


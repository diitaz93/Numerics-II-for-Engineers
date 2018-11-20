function [ xh, Lh, fh ] = a03ex03getBVP( p )
% Assignment 03, Exercise 03b, by Shane Alpert (404579), 
% Juan Sebastian Diaz (405385), Bo Yan (403787)
% 
% Sets up the grid xh, the sparse matrix Lh, and the right hand side fh of 
% the corresponding linear system for the refinement level n=2^p-1. This 
% discretized problem is later solved (and the error determined) in the 
% % function "a03ex03solve".
N = 2^p - 1;    % defines number of internal grid points
h = 1/(N+1);    % defines mesh size/spacing
xh = (linspace(0,1,N+2));
xh = xh(2:end-1)';
% defines grid points
% Creation of Lh sparse matrix
% [ (diagonal) (upper off-diagonal) (lower off-diagonal)]
r = [1:(N) 1:(N-1) 2:(N)];
% [ (diagonal) (upper off-diagonal) (lower off-diagonal)]
c = [1:(N) 2:(N) 1:(N-1)];
% [ (diagonal) (upper off-diagonal) (lower off-diagonal)]
v = [((2/(h^2))+1)*ones(1,N) -((1+2*h)/(h^2))*ones(1,N-1) ((2*h-1)/(h^2))*ones(1,N-1)];
Lh = sparse (r,c,v,N,N);
% Inhomogeneity
fh = -3*xh.^3 + 40*xh.^2 - 14*xh - 7;
fh(1)=fh(1)+(1-2*h)/h^2;
fh(N)=fh(N)+(2+4*h)/h^2;
end


function [t] = a02e04getPDEType(A, b)
% Assignment 02, Exercise 04, by Shane Alpert, Juan Sebastian Diaz, Bo Yan
%
% A02E04GETPDETYPE returns the type of a constant coefficient PDE given the
% second order coefficient matrix A and the first order coefficient vector
% b. t takes the value 1 for eliptic PDEs, 2 for parabolic PDEs, 3 for
% hyperbolic PDEs, and 4 for unclassified PDEs.
%   The function uses the eigenvalues' signs and the rank of [AB] to
%   classify the PDE system.
tic
ev = eig(A); % Eigenvalues of A
T = [A,b]; % Concatenation of matrix A and vector b
s = sign(ev); % Vector of signs of eigenvalues
n = length(A(:,1)); % Dimension of matrix A
% If statement for elliptic 
% No non-zero elements and all of the same sign
if (nnz(s)==n && length(unique(s))==1) 
    t = 1;
% If statement for hyperbolic
% No non-zero elements and n-1 eigenvalues of same sign
elseif (nnz(s)==n && (sum(s==1)==(n-1) || sum(s==-1)==(n-1) ))
        t = 3;
% If statement for parabolic
% n-1 non-zero elements, rank of [AB] diferent to zero and all non zero
% eigenvalues of the same sign
elseif (rank(T)==n && nnz(s)==(n-1) && length(unique(s))==2)
    t = 2;
else
    t = 0;
end
toc
    
end


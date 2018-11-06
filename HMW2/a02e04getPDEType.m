function [t] = a02e04getPDEType(A, b)
% Assignment 02, Exercise 04, by Shane Alpert, Juan Sebastian Diaz, Bo Yan
%
% A02E04GETPDETYPE returns the type of a constant coefficient PDE given the
% second order coefficient matrix A and the first order coefficient vector
% b. t takes the value 1 for eliptic PDEs, 2 for parabolic PDEs, 3 for
% hyperbolic PDEs, and 4 for unclassified PDEs.
%   The function uses the eigenvalues...
tic
ev = eig(A);
T = [A,b];
s = sign(ev);
n = size(A(:,1));
if (nnz(s)==0 && length(unique(s))==1)
    t = 1;
elseif (nnz(s)==0 && (sum(s==1)==(n-1) || sum(s==-1)==(n-1)) )
        t = 3;
elseif (rank(T)~=0 && nnz(s)==1 && length(unique(s))==2)
    t = 2;
else
    t = 0;
end
toc
    
end


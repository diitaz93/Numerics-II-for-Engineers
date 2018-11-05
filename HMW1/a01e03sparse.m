function Kn = a01e03sparse(n)
% Assignment 01, Exercise 03, by Shane Alpert, Juan Sebastian Diaz, Bo Yan
%
% A01E03SPARSE(n) returns the sparse matrix Kn of dimension n with entries
% 2 in the diagonal, -1 in the nearest off-diagonals and zero elswhere.

% Create array with row indexes of non-zero entries of matrix
% [ (diagonal) (upper off-diagonal) (lower off-diagonal)] 
r = [1:n 1:n-1 2:n];
% Create array with column indexes of non-zero entries of matrix
% [ (diagonal) (upper off-diagonal) (lower off-diagonal)]
c = [1:n 2:n 1:n-1];
% Create array with values of non-zero entries of matrix
% [ (diagonal) (off-diagonal (upper & lower))]
v = [2*ones(1,n) -1*ones(1,2*n-2)];
% Use sparse function with created arrays as parameters, also allocate
% memory for a size n-by-n matrix
Kn = sparse(r,c,v,n,n);
end

% What are the advantages of working with sparse matrices?
% For big matrices, using sparse matrices improves efficiency of memory
% usage due to the fact that memory is not wasted saving huge amounts of
% zeros, but only the non-zero values. Sparse matrices also save time when
% doing matrix operations such as mulitplication,
% etc.


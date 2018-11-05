% Assignment 01, Exercise 04 a), by Shane Alpert, Juan Sebastian Diaz,
% Bo Yan
%% Exercise 4a) Function
function x = a01e04invKn( b )

%A01E04INVKN computes the solution x of Kn*x=b for any given b vector,
% being Kn the matrix used in the previous exercises.
% The function uses the Cholesky decomposition of matrix Kn, doing forward
% substitution for the lower triangular matrix (C') and backward
% substitution for the upper triangular matric (C). The values of C and C'
% were obtained by a derivation of the general form of the Cholesky
% decomposition for matrices of the type Kn. 

n = length(b); % Dimension of vector
% First, we solve the problem C'*y=b (solve for y)
y = ones(n,1); % Declaration of y vector
y(1) = (b(1)/sqrt(2)); % Initialization of first entry

% Forward substitution
for j=2:n
    % General form of forward substitution using Cholesky decomposition
    % for Kn matrix.
    y(j) = b(j)*sqrt(j/(j+1)) + y(j-1)*sqrt((j-1)/(j+1));
end
% Second, we solve the problem C*x=y (solve for x)
x = ones(n,1);% Declaration of y vector
x(n) = y(n) * sqrt(n/(n+1));% Initialization of last entry

% Backward substitution (C*x=y)
for j=n-1:-1:1
    % General form of backward substitution using Cholesky decomposition
    % for Kn matrix.
    x(j) = sqrt(j/(j+1))*y(j) + j*x(j+1)/(j+1);
end
end
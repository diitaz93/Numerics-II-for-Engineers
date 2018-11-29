function [Lh,fh] = a05e04getPDE(p,i)
% a05e04getPDE Sets up sparse matrix Lh and the right hand side fh of the 
% linear system Lh*uh=uh for the refinement level p on the domain
% (0,1) * (0,1). We use the standard five-point stencil on a uniform mesh
% with lexicographical order. 

% Initialize grid spacing and mesh size
N = (2^p) - 1;
h = 1/(N+1);
x = (1:N)*h;
[X,Y] = meshgrid(x);

% Creates fh
if i == 1
    fh = 6*X.*Y.*(2-Y.^2-X.^2);
elseif i == 2
    fh = 10*(pi)^2*sin(3*pi*X).*sin(pi*Y);
else
    msgbox('Invalid Value: Please enter 1 or 2.', 'Error','error');    
end

% Transposes fh and turns it into a vector using lexicographical ordering
fh = fh';
fh = fh(:);

% Set-up sparse matrix S for use in creating Lh
% Create array with row indices of non-zero entries of matrix
% [ (diagonal) (upper off-diagonal) (lower off-diagonal)] 
r = [1:N 1:N-1 2:N];
% Create array with column indices of non-zero entries of matrix
% [ (diagonal) (upper off-diagonal) (lower off-diagonal)]
c = [1:N 2:N 1:N-1];
% Create array with values of non-zero entries of matrix
% [ (diagonal) (off-diagonal (upper & lower))]
v = [2*ones(1,N) -1*ones(1,2*N-2)];
% Use sparse function with created arrays as parameters, also allocate
% memory for a size n-by-n matrix
S = sparse(r,c,v,N,N);

% Define Lh using Kronecker products
Lh = - ( kron(eye(N),S) + kron(S,eye(N)) ) / (h^2);


end


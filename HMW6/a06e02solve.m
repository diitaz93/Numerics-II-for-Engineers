function [lambda,error,eoc] =  a06e02solve ()
% a06e02solve calulates solves the numerically the PDE given in the
% excercise with function handles defined below. It returns the coefficient
% lambda of the extended system, the error between the numeric and exact
% solution, and the experimental order of convergence.

% Assignment 06, Exercise 03, by Shane Alpert (404579), 
% Juan Sebastian Diaz (405385), Bo Yan (403787) 

%% Initialization
% Define the function handles given in excercise for f, g and u
f = @(x,y)cos(2*pi*x).*exp(y.^3).*(4*pi^2-6*y-9*y.^4);
g = @(x,y) 3*cos(2*pi*x).*exp(y).*(y==1);
u = @(x,y) cos(2*pi*x).*exp(y.^3);

lambda = zeros(8,1);
errors = zeros(8,1);
h = zeros(8,1);
eoc = zeros(7,1);
%% Loop over the values of p
for p=1:8
    % Determination of grid parameters
    N=2^p-1;
    h(p)=1/(N+1);
    %% Lh,fh
    if p==1
        Lh = 0;
        X = 1/2;Y=1/2;

        % Column 1, all rows, i=1
        x1 = 0;
        y1 = 1/2;

        g1 = g(x1,y1);

        % Row 1, all columns,j=1
        x2 = 1/2;
        y2 = 0;

        g2 = g(x2,y2);

        % Column N, all rows, i=N
        x3 = 1;
        y3 = 1/2;

        g3 = g(x3,y3);

        % Row N, all columns, j=N
        x4 = 1/2;
        y4 = 1;

        g4 = g(x4,y4);

        Rhf= f(X,Y);
        G=(g1 + g2 + g3 + g4)/h(p);
        fh = Rhf+ G;
        fh = fh';
        fh = fh(:);
    else
    [Lh,fh] = a06e02getPDE(N,f,g);
    end
    %% uh
    Lht=[Lh ones(N^2,1) ; ones(1,N^2) 0];
    fht = [fh ; 0];
    uht = Lht\fht;
    uh = uht(1:end-1);
    %% Lambda
    lambda(p) = uht(end);
    %% Exact solution
    [X,Y] = meshgrid(h(p):h(p):1-h(p));
    Rhu = u(X,Y)';
    Rhu = Rhu(:);
    %% Error
    errors(p) = norm( (Rhu-uh) , Inf);
    %% EOC 
    if p>1
    eoc(p-1) = ( log(errors(p)) - log(errors(p-1)) ) / ...
                ( log(h(p)) - log(h(p-1)) ) ;
    end
end
%% Monitor EOC
plot((2:8)',eoc);
xlabel('p')
title('Experimental order of convergence')
end


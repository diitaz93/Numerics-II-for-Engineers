function [errors] = a05e04solve(i)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Initialize p vector
p=1:9;
errors = zeros(1,9);
h = zeros(1,9);
EOC = zeros(1,8);

if i == 1
    
    for j = 1:9
        % Initialize grid spacing and mesh size
        N = (2^p(j)) - 1;
        h(j) = 1/(N+1);
        x = (1:N)*h(j);
        [X,Y] = meshgrid(x);
        
        % Exact solution
        u = X.*Y + (X.^3).*(Y.^3) - (X.^3).*(Y) - (X).*(Y.^3);
        u = u';
        u = u(:);
        
        % Numerical approximation
        [Lh,fh] = a05e04getPDE(p(j),i);
        uh = Lh\fh;
        
        % Determine error between exact solution and approximation in max
        % norm
        errors(j) = norm( (u-uh) , Inf); 
        
        
        % Determine the experimental order of convergence
        if j==1
            continue;
        else
            EOC(j-1) = ( log(errors(j)) - log(errors(j-1)) ) / ...
                ( log(h(j)) - log(h(j-1)) ) ;
        end

    end
    
    loglog(h,errors)
    table((p(2:end))',(EOC)','VariableNames',{'p','EOC'})
    
elseif i == 2
    
    for j = 1:9
        % Initialize grid spacing and mesh size
        N = (2^p(j)) - 1;
        h(j) = 1/(N+1);
        x = (1:N)*h(j);
        [X,Y] = meshgrid(x);
        
        % Exact solution
        u = sin(3*pi*X).*sin(pi*Y);
        u = u';
        u = u(:);
        
        % Numerical approximation
        [Lh,fh] = a05e04getPDE(p(j),i);
        uh = Lh\fh;
        
        % Determine error between exact solution and approximation in max
        % norm
        errors(j) = norm( (u-uh) , Inf); 
        
        % Determine the experimental order of convergence
        if j==1
            continue;
        else
            EOC(j-1) = ( log(errors(j)) - log(errors(j-1)) ) / ...
                ( log(h(j)) - log(h(j-1)) ) ;
        end

    end
    
    loglog(h,errors)
    table((p(2:end))',(EOC)','VariableNames',{'p','EOC'})

else 
    msgbox('Invalid Value: Please enter 1 or 2.', 'Error','error');
    

end
   

end


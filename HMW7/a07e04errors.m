u = @(x) x.^4-2*x.^2+0.5;
errors = ones(1,8);
for p=2:9
    N = 2^p-1;
    h = 2/(N+1);
    xh = [-1+h:h:1-h]';
    Rhu = u(xh);
    [Lh, fh] = a07e04getPDE(p,0);
    uh = Lh\fh;
    errors(p-1)=norm( (Rhu-uh) , Inf);
end
loglog(2:9,errors)
title('Infinity norm of difference between exact solution and computation')
ylabel('Error')
xlabel('p')
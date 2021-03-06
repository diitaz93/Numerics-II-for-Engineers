p = 9;
beta = [-1:0.4:1 0];
N = 2^p-1;
h = 2/(N+1);
x = -1+h:h:1-h;
figure()
for i =1:length(beta)
    [Lh,fh] = a07e04getPDE(9,beta(i));
    uh = Lh\fh;
    plot(x,uh);
    hold on
end
title('$u_h$ for different Betas','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$u_h$','Interpreter','latex')
legend(['\beta =' num2str(beta(1))],...
       ['\beta =' num2str(beta(2))],...
       ['\beta =' num2str(beta(3))],...
       ['\beta =' num2str(beta(4))],...
       ['\beta =' num2str(beta(5))],...
       ['\beta =' num2str(beta(6))],...
       ['\beta =' num2str(beta(7))],...
       'location', 'best')
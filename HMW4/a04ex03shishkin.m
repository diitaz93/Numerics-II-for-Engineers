function xh=a04ex03shishkin(N,sigma)
H = (1-sigma)/N;
h = sigma/N;
i1 = [2:N]';
i2=[N+1:2*N]';
xh =[(i1-1)*H; (1-sigma) + (i2-N-1)*h];
% xh=[linspace(0,(N-1)*(1-sigma)/N,N) linspace(1-sigma,1,N+1)];
% xh=xh'
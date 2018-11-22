function uh=a04ex03solve(eps,xh,flag)
%function uh=a04ex03solve(eps,xh,flag)
%xh:x1,...,xn
%flag:-,+,0
consts=[eps,1,0,0,0];
f=@(x) 0.*x+1;
[Lh,fh]=a04ex02getPDE(xh,f,consts,flag);
uh=Lh\fh;
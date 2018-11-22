function [Lh,fh]=a04ex02getPDE(xh,f,consts,flag)
% A04EX02GETPDE returns the Lh matrix and the fh vector of a second order
% PDE system with boundary conditions. 
%xh:x1,...,xn
%consts=[a,b,c,al,be];-au"+bu'+cu=f;al at  0;be at 1
%flag -,+,0
a=consts(1);
b=consts(2);
c=consts(3);
al=consts(4);
be=consts(5);

n=size(xh,1);
%fh=fh1+fh2;f(xh)->fh1;a,b,c,alpha,beta->fh2;
fh1=f(xh);

%->h:(n+1)*1
xh=[xh;1];% xh->x1,...,xn,beta
temp=[0;xh(1:end-1)];%temp->alpha,x1...,xn
h=xh-temp;%h->h1,...,h(n+1);hi=xi-x(i-1)
%L=-aA+bB+cI
     %A:
     Da=-2./h(1:end-1);Da=Da./h(2:end);
     Ua=2./h(2:end-1);Ua=Ua./(h(1:end-2)+h(2:end-1));
     La=2./h(2:end-1);La=La./(h(2:end-1)+h(3:end));
     Va=[Da;Ua;La];
     CL=[1:n 2:n 1:(n-1)]; CL=CL';
     RL=[1:n 1:(n-1) 2:n];RL=RL';
     A=sparse(RL,CL,Va,n,n);
     %B:
switch flag
    case '-'
     Db=1./h(1:end-1);
     Lb=-1./h(1:end-2);
     Vb=[Db;zeros(n-1,1);Lb]; 
     
     fh2=[2*a*al/h(1)/(h(1)+h(2))+b*al/h(1);zeros(n-2,1);2*a*be/h(end)/(h(end-1)+h(end))];
    case '+'
     Db=-1./h(2:end);
     Ub=1./h(2:end-1);
     Vb=[Db;Ub;zeros(n-1,1)];
     
     fh2=[2*a*al/h(1)/(h(1)+h(2));zeros(n-2,1);2*a*be/h(end)/(h(end-1)+h(end))-b*be/h(end)];
    case '0'
     Ub=1./(h(1:end-2)+h(2:end-1));
     Lb=1./(h(2:end-1)+h(3:end));
     Vb=[zeros(n,1);Ub;Lb];
     
     fh2=[2*a*al/h(1)/(h(1)+h(2))+b*al/(h(1)+h(2));zeros(n-2,1);2*a*be/h(end)/(h(end-1)+h(end))-b*be/(h(end-1)+h(end))];
    otherwise
        disp('Unknown flag');
end
     B=sparse(RL,CL,Vb,n,n);
Lh=-a.*A+b.*B+c.*eye(n,n);
fh=fh1+fh2;
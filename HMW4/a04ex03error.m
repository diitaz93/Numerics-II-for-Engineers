function [err,uex]=a04ex03error(eps,xh,uh)
f=@(x) x-((exp(-(1-x)/(eps)))-(exp(-1/eps)))/(1-exp(-1/eps));
uex=f(xh);
err=norm(uh-uex,Inf);
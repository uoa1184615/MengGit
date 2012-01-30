function dhudt=hol_matrix(t,hu)
global n m b gam dx ph pu qh qu qhl qhr qul qur qhc quc l nu1 nu2

uu=nan(n+1,m+2);
uu=uu(:);
uu(qh)=hu(ph);
uu(qu)=hu(pu);

%periodic boundary condition.
uu(b:n)=uu(end-l-4:end-l-1);
uu(end-l+1:end-b)=uu(n+2:n+b+1);

% insert couple boundary conditions on each tooth by 
% applying periodic bcs
uu(qhr)=(1+gam)/2*uu(qhc)+(1-gam)/2*uu(quc-l);
uu(qhl)=(1-gam)/2*uu(qhc)+(1+gam)/2*uu(quc-l);
uu(qur)=(1+gam)/2*uu(qhc+l)+(1-gam)/2*uu(quc);
uu(qul)=(1-gam)/2*uu(qhc+l)+(1+gam)/2*uu(quc);
uu(n)=uu(end-l-b);
uu(end-n)=uu(l+4);

% insert values for the second derivative
%uu(quc-b)=0;
%uu(quc+b)=0;

dhudt=nan(length(hu),1);
dhudt(ph)=-(uu(qh+1)-uu(qh-1))/dx/2;
dhudt(pu)=-(uu(qu+1)-uu(qu-1))/dx/2-nu1*uu(qu);
%dhudt(pu)=-(uu(qu+1)-uu(qu-1))/dx/2-nu1*uu(qu) ...
               % +nu2*(uu(qu+2)-2*uu(qu)+uu(qu-2))/dx^2;





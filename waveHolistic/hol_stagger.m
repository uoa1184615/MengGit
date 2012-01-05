function dhudt=hol_stagger(t,hu)
global n m b gam dx ph pu qh qu qhl qhr qul qur qhc quc l nu

uu=nan((n+1)*(m+1)+1,1);
uu(qh)=hu(ph);
uu(qu)=hu(pu);

%periodic boundary condition.
uu(1:b)=uu(end-l:end-b-1);
uu(end-b+1:end)=uu(b+2:n+2);

% insert couple boundary conditions on each tooth by 
% applying periodic bcs
uu(qhr)=(1+gam)/2*uu(qhc)+(1-gam)/2*uu(quc-l);
uu(qhl)=(1-gam)/2*uu(qhc)+(1+gam)/2*uu(quc-l);
uu(qur)=(1+gam)/2*uu(qhc+l)+(1-gam)/2*uu(quc);
uu(qul)=(1-gam)/2*uu(qhc+l)+(1+gam)/2*uu(quc);

% insert values for the second derivative
uu(quc-b)=0;
uu(quc+b)=0;

dhudt=nan(length(hu),1);
dhudt(ph)=-(uu(qh+1)-uu(qh-1))/dx/2;
dhudt(pu)=-(uu(qu+1)-uu(qu-1))/dx/2-nu*uu(qu);
%dhudt(pu)=-(uu(qu+1)-uu(qu-1))/dx/2-nu*uu(qu) ...
%                +nu*(uu(qu+2)-2*uu(qu)+uu(qu-2))/dx^2;





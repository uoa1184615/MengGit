global n m b gam dx ph pu qh qu qhl qhr qul qur qhc quc l
n=19; %n=7,11,...
m=12; % m=even
b=(n+1)/2;
l=n+1;
L=2*pi;
D=L/m;
d=D/2;
gam=d/D/2;
dx=d/(n-1);

% define indexes into the dynamic variables
p=1:m*(n-2);
ph=p(2:2:end)';
pu=p(1:2:end)';

% q used for equations, shifted for BCs
q=zeros(n+1,m); q(:)=1+b+(1:(n+1)*m);
qh1=q(3:2:end-2,1:2:m); qh1=qh1(:);
qh2=q(2:2:end-2,2:2:m); qh2=qh2(:);
qh=[qh1;qh2]; qh=sort(qh);
qu1=q(2:2:end-2,1:2:m); qu1=qu1(:);
qu2=q(3:2:end-2,2:2:m); qu2=qu2(:);
qu=[qu1;qu2]; qu=sort(qu);
qhl=q(1,1:2:m)'; qhr=q(end-1,1:2:m)';
qul=q(1,2:2:m)'; qur=q(end-1,2:2:m)';
quc=q(b,1:2:m)'; qhc=q(b,2:2:m)';

% define the x-axis
x=nan((m+1)*(n+1)+1,1);
x(1:b)=(0:b-1)*dx;
for i=1:m
    x((n+1)*i+(-b+2:b))=D*i+(-b+1:b-1)*dx;
end
x(end-b+1:end)=L+D+(-b+1:0)*dx;

% initial values
hu=nan(m*(n-2),1);
hu(ph)=1+0.1*sin(x(qh));
hu(pu)=0*0.1*sin(x(qu)+pi/2);

%ts0=linspace(0,10,4);
[ts hus]=ode15s(@hol_stagger,[0 10],hu);

% plot the results
hs=nan(length(ts),(m+1)*(n+1)+1);
hs(:,qh)=hus(:,ph);
%hs(:,1:2:b)=hus(:,end-b+1:2:end);
%hs(:,end-b+1:2:end)=hus(:,1:2:b);
%hs(:,qhr)=(1+gam)/2*hs(:,qhc)+(1-gam)/2*hs(:,quc-l);
%hs(:,qhl)=(1-gam)/2*hs(:,qhc)+(1+gam)/2*hs(:,quc-l);
us=nan(length(ts),(m+1)*(n+1)+1);
us(:,qu)=hus(:,pu);
%us(qur)=(1+gam)/2*us(qhc+l)+(1-gam)/2*us(quc);
%us(qul)=(1-gam)/2*us(qhc+l)+(1+gam)/2*us(quc);
return
subplot(2,1,1);plot(x,hs(end,:),'ko');
clf();
subplot(2,1,1);plot(x,hs,'o');xlabel('x');ylabel('h');
subplot(2,1,2);plot(x,us,'o');xlabel('x');ylabel('u');


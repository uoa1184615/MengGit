global n m b gam dx ph pu qh qu qhl qhr qul qur qhc quc l nu1 nu2
n=7; %n=7,11,...
m=4; % m=even
b=(n+1)/2;
l=n+1;
L=10;
D=L/m;
d=D/2;
gam=d/D/2;
nu1=0.1;
nu2=1;
dx=d/(n-1);

% define indexes into the dynamic variables
p=1:m*(n-2);
ph=p(2:2:end)';
pu=p(1:2:end)';

% q used for equations, shifted for BCs
q=zeros(n+1,m+2); q(:)=1:(n+1)*(m+2);
qh1=q(3:2:end-2,2:2:end-1); qh1=qh1(:);
qh2=q(2:2:end-2,3:2:end-1); qh2=qh2(:);
qh=[qh1;qh2]; qh=sort(qh);
qu1=q(2:2:end-2,2:2:end-1); qu1=qu1(:);
qu2=q(3:2:end-2,3:2:end-1); qu2=qu2(:);
qu=[qu1;qu2]; qu=sort(qu);
qhl=q(1,2:2:end-1)'; qhr=q(end-1,2:2:end-1)';
qul=q(1,3:2:end-1)'; qur=q(end-1,3:2:end-1)';
quc=q(b,2:2:end-1)'; qhc=q(b,3:2:end-1)';

% define the x-axis
x=nan(n+1,m+2);
x=x(:);
x(1:b-1)=0;
x(b:n)=(0:b-1)*dx;
for i=1:m+1
    x((n+1)*i+(1:n))=D*i+(-b+1:b-1)*dx;
end
%x(end-b+1:end)=L+D+(-b+1:0)*dx;

% initial values
hu=nan(m*(n-2),1);
hu(ph)=1+0.1*sin(2*pi/L*x(qh));
hu(pu)=0.1*sin(2*pi/L*x(qu)+pi/2);

ts0=linspace(0,40,100);
[ts hus]=ode15s(@hol_matrix,ts0,hu);

% plot the results
for it=1:length(ts)
hs=nan(length(ts),(m+2)*(n+1));
hs(it,qh)=hus(it,ph);
us=nan(length(ts),(m+2)*(n+1));
us(it,qu)=hus(it,pu);
figure(1);
    subplot(2,1,1),plot(x,hs,'bo');
       xlabel('x');ylabel('h');
       axis([0 L+D 0.8 1.2]);
       title(['t=' num2str(ts(it))]);
    subplot(2,1,2),plot(x,us,'bo');
       xlabel('x');ylabel('u');
       axis([0 L+D -0.2 0.2]);
       title(['t=' num2str(ts(it))]);
       drawnow
       %pause
       M(:,it)=getframe(gcf);
end
%movie(M);
map=colormap; 
mpgwrite(M, map,'filename');


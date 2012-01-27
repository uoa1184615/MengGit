global n m b gam dx ph pu qh qu qhl qhr qul qur qhc quc l nu
n=7; %n=3,7,11,...
m=10; % m=even
b=(n+1)/2;
l=n+1;
L=2*pi;
D=L/m;
d=D/4;
gam=d/D/2;
nu=0.01;
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
k=2*pi/L;
hu=nan(m*(n-2),1);
hu(ph)=1+0.1*sin(k*x(qh));
hu(pu)=  0.1*sin(k*x(qu));

ts0=linspace(0,20,100);
if not(exist('OCTAVE_VERSION','builtin'))
[ts hus]=ode15s(@hol_stagger,ts0,hu);
end;

% plot the results
for it=1:length(ts)
hs=nan(1,(m+1)*(n+1)+1);
hs(qh)=hus(it,ph);
us=nan(1,(m+1)*(n+1)+1);
us(qu)=hus(it,pu);
figure(1);
    subplot(2,1,1),plot(x,hs,'bo');
       xlabel('x');ylabel('h');
       axis([0 L+D 0.8 1.2]);
    subplot(2,1,2),plot(x,us,'bo');
       xlabel('x');ylabel('u');
       axis([0 L+D -0.1 0.1]);
       drawnow
       pause(0.1)
       M(:,it)=getframe(gcf);
end
movie(M);
map=colormap; 
mpgwrite(M, map,'filename');


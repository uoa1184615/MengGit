global n m b gam dx ph pu qh qu qhl qhr qul qur qhc quc l nu
% from these eigenvalues it certainly looks like there
% should be travelling waves
n=7; %n=3,7,11,...
m=6; % m=even, but mod(m,4)=2 may be best
b=(n+1)/2;
l=n+1;
L=2*pi;
D=L/m;
d=D/4; % the size of the patches
gam=d/D/2;
nu=0.01;
dx=d/(n-1);
eps=1e-5;

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

% check the equilibrium
hu=nan(m*(n-2),1);
hu(ph)=1;
hu(pu)=0;
equil=norm(hol_stagger(0,hu))

N=length(pu)+length(ph);
J=nan(N,N);
for jj=1:N
    J(:,jj)=hol_stagger(0,hu+eps*(jj==(1:N))');
end
J=J/eps;
[v,d]=eig(J);
d=diag(d);
% order smallest to largest in magnitude
[dd,j]=sort(abs(d)); 
d=d(j), v=v(:,j);
% plot eigenvalues
figure(1), clf()
%horiz=1*asinh(real(d));
%vert=1*asinh(imag(d));
horiz=asinh(real(d))+1e-4*(2*rand(size(asinh(real(d))))-1);
vert=asinh(imag(d))+1e-4*(2*rand(size(asinh(imag(d))))-1);
plot(horiz,vert,'o');xlabel('asinh(Re)');ylabel('asinh(Im)');
% plot eigenvectors
hus=v';
maxv=max(abs(v(:)));
for it=2:2:min(10,size(hus,2))
    hs=nan(1,(m+1)*(n+1)+1);
    hs(qh)=hus(it,ph);
    us=nan(1,(m+1)*(n+1)+1);
    us(qu)=hus(it,pu);
    figure(it),clf()
    subplot(2,1,1),plot(x,real(1*hs),'bo',x,imag(1*hs),'r+');
       xlabel('x');ylabel('h');
       axis([0 L+D -maxv maxv]);
    subplot(2,1,2),plot(x,real(1*us),'bo',x,imag(1*us),'r+');
       xlabel('x');ylabel('u');
       axis([0 L+D -maxv maxv]);
       pause(1)
end

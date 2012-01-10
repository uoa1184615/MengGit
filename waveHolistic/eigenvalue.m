global n m b gam dx ph pu qh qu qhl qhr qul qur qhc quc l nu
n=7; %n=7,11,...
m=4; % m=even
b=(n+1)/2;
l=n+1;
L=2*pi;
D=L/m;
d=D/2;
gam=d/D/2;
nu=0;
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
d=diag(d)
horiz=10e16*asinh(real(d));
vert=1*asinh(imag(d));
%horiz=asinh(real(d))+0.1*(2*rand(size(asinh(real(d))))-1);
%vert=asinh(imag(d))+0.1*(2*rand(size(asinh(imag(d))))-1);
plot(horiz,vert,'o');%xlabel('asinh(Re)');ylabel('asinh(Im)');
global n m r dx nu0 nu2 gam
n=9;  % 5,9,... number of microscale grid points, not incl edges
m=6; % even number of patches
r=1/6 ;% ratio of patch size to macrogrid size
L=2*pi; % length of domain
D=L/m;
d=2*r*D;
dx=d/(n+1); % patch size includes edges, to give microgrid size
nu0=0.1;
nu2=0.03;
gam=1;
eps=1e-5; % order parameter

% define the x-coords, omit patch edges
[hu,x]=meshgrid((0:m-1)*D,(0:n-1)*dx);
x=x+hu;

% check the equilibrium
hu(1:2:end)=1;
hu(2:2:end)=0;
equil=norm(hol_staggermatrix(0,hu))

J=nan(n*m,n*m);
for jj=1:n*m
    J(:,jj)=hol_staggermatrix(0,reshape(hu(:)+eps*(jj==(1:n*m))',n,m));
end

J=J/eps;
[v,d]=eig(J);
d=diag(d);
% order smallest to largest in magnitude
[dd,j]=sort(abs(d)); 
d=d(j), v=v(:,j);
% plot eigenvalues
figure(1), clf()
horiz=1*asinh(real(d));
vert=1*asinh(imag(d));
%horiz=asinh(real(d))+1e-4*(2*rand(size(asinh(real(d))))-1);
%vert=asinh(imag(d))+1e-4*(2*rand(size(asinh(imag(d))))-1);
%plot(horiz,vert,'o');xlabel('asinh(Re)');ylabel('asinh(Im)');
plot(horiz,imag(d),'o');xlabel('asinh(Re)');ylabel('Im');
% plot eigenvectors
return

hus=v';
maxv=max(abs(v(:)));
for it=2:2:min(10,size(hus,2))
    hs=nan(1,n*m);
    hs(1,1:2:end)=hus(it,1:2:end);
    us=nan(1,n*m);
    us(1,2:2:end)=hus(it,2:2:end);
    figure(it),clf()
    subplot(2,1,1),plot(x(:),real(1*hs),'bo',x(:),imag(1*hs),'r+');
       xlabel('x');ylabel('h');
       axis([0 L+D -maxv maxv]);
    subplot(2,1,2),plot(x(:),real(1*us),'bo',x(:),imag(1*us),'r+');
       xlabel('x');ylabel('u');
       axis([0 L+D -maxv maxv]);
       pause(1)
end

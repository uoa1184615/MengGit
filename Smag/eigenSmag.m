global k kp0 km0 k0p k0m dx dy theta g b c
nx=20;
ny=45;    
g=1;
Lx=30;
Ly=20;
dx=Lx/nx;
dy=Ly/ny;
c=0.1;
theta=0.003;
mag=0.7;
alpha=9;
beta=3;
i=1:nx;
j=1:ny;
x=(i-1)*dx; 
y=(j-1)*dy;
[y,x]=meshgrid(y,x);

k=zeros(nx,ny);
k(:)=1:nx*ny;
kkk=[k k k;k k k;k k k];
kp0=kkk(nx+i+1,ny+j);
km0=kkk(nx+i-1,ny+j);
k0p=kkk(nx+i,ny+j+1);
k0m=kkk(nx+i,ny+j-1);
k=k(:);
kp0=kp0(:);
km0=km0(:);
k0p=k0p(:);
k0m=k0m(:);
% description of curved bed
A=0.7;
f=A*cos(2*pi/Lx*x);
b1=1-max(0,1-(y-f-alpha).^2/beta^2).^2;
 b=mag*b1;

% check the equilibrium
h=ones(nx,ny)-b;
u=sqrt(0.985*sin(theta).*h/0.0029);
v=zeros(nx,ny);
huv=[h(:)';u(:)';v(:)'];
huv=huv(:);
equil=norm(smagcurthree(0,huv))

return
J=nan(3*nx*ny,3*nx*ny);
for jj=1:3*nx*ny
    J(:,jj)=smagcurthree(0,huv(:)+eps*(jj==(1:3*nx*ny))');
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
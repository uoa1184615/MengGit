global k kp0 km0 k0p k0m dx dy theta g b c
nx=20;
ny=45;    
g=1;
Lx=30;
Ly=20;
dx=Lx/nx;
dy=Ly/ny;
c=0.1;
theta=0.005;
mag=0.8;
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
A=0.8;
f=A*cos(3*pi/Lx*x);
b1=1-max(0,1-(y-f-alpha).^2/beta^2).^2;
 b=mag*b1;
%b=mag*ones(nx,ny);
 %b(:,3*ny/10:7*ny/10)=0;
 %b(:,2:5*ny/10)=0;
 %%%%curved channel
return
 
 %initial values
h=0*0.01*sin(2*pi/Lx*x).*sin(2*pi/Ly*y)-b+1;
u=0*sqrt(0.985*sin(theta).*h/0.0029) ...
    +0*0.001*rand(nx,ny)+0.001*sin(2*pi/Lx*x);
v=zeros(nx,ny);
huv=[h(:)';u(:)';v(:)'];
huv=huv(:);
%return

%ts0=linspace(0,30,4);
%[ts huvs]=ode15s(@smagcurthree,ts0,huv);
[ts huvs]=ode15s(@smagcurthree,[0 100],huv);
%return
clf()
for it=1:length(ts)
    huv=huvs(it,:);
    h(:)=huv(1:3:end);
    u(:)=huv(2:3:end);
    v(:)=huv(3:3:end);
    figure(1);
    subplot(2,2,1),surf(x,y,h+b);
    %plot(x,h(:,end/2)+B(:,end/2),'g');
       xlabel('x');ylabel('y');zlabel('free surface');
       axis([0 Lx 0 Ly 0.2 1.8]);
        title(['t=' num2str(ts(it))]);
    subplot(2,2,2),surf(x,y,b);
       xlabel('x');ylabel('y');zlabel('bed surface');
       %axis([0 Lx 0 Ly 0 1]);
    subplot(2,2,3),surf(x,y,u);
       xlabel('x');ylabel('y');zlabel('mean u');
        axis([0 Lx 0 Ly 0 3]);
         title(['t=' num2str(ts(it))]);
        %axis([0 Lx 0 Ly -0.5 1]);
    subplot(2,2,4),surf(x,y,v);
        xlabel('x');ylabel('y');zlabel('mean v');
       axis([0 Lx 0 Ly -0.5 0.5]);
        title(['t=' num2str(ts(it))]);
        %axis([0 Lx 0 Ly -0.5 1]);            
 drawnow
   %pause;
     M(:,it)=getframe(gcf);
end
%movie(M);
map=colormap; 
mpgwrite(M, map,'filename');
%qq=sqrt(u.^2+v.^2);
%qu=[u(:,end),u,u(:,1)];
%qu=[qu(end,:);qu;qu(1,:)];
%qv=[v(:,end),v,v(:,1)];
%qv=[qv(end,:);qv;qv(1,:)];
%fs=h+B;
%qh=[fs(:,end),fs,fs(:,1)];
%qh=[qh(end,:);qh;qh(1,:)];
%qi=2:nx+1;
%qj=2:ny+1;
%txz=0.002784*u.*qq+0.008615*h*g*sin(theta) ...
 %   -0.008615*h*g*cos(theta)/2.*(qh(qi+1,2:end-1)-qh(qi-1,2:end-1));
%tyz=0.002784*v.*qq ...
 %   -0.008382*h*g*cos(theta)/2.*(qh(2:end-1,qj+1)-qh(2:end-1,qj-1));
%txzf=-0.0001467*u.*qq-0.006461*h*g*sin(theta) ...
    %+0.006461*h*g*cos(theta)/2.*(qh(qi+1,2:end-1)-qh(qi-1,2:end-1));


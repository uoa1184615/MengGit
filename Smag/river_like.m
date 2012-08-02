% Simulations in river-like channel.
% MC and AJR, 19 July 2012

global k kp0 km0 k0p k0m dx dy theta g b gam
nx=40; % grids in x-direction 
ny=40; % grids in y-direction   
g=1; % gravity
gam=1; % coupling parameter
Lx=40; % considered length in x-direction
Ly=20; % considered length in y-direction
dx=Lx/nx; % spatial step in x-direction
dy=Ly/ny; % spatial step in y-direction
theta=0.001; % slope in x-direction

% define the coordinate and the river-like channel.
% multi-D linear spline prediction
xk=[.2 .2 .2 .2 .5 .5 .5 .5 .25 .25 .25 .25 .4 .4 .4 .4 .7 .7 .7 .7 .8 .8 .8 .8 .9 .9 .9 .9];
yk=[.1 .3 .6 .9 .3 .5 .7 .9 .1 .4 .6 .9 .3 .4 .5 .9 .1 .4 .6 .9 .1 .3 .45 .7 .1 .3 .6 .9];
fk=[.8  0  .8  .8  .8  .8  0  .8 .8 0 .8 .8 .8 .8 0 .8 .8 .8 0 .8 .8 .8 0 .8 .8 0 .8 .8];
% force column vector and make periodic??
n=length(xk);
xk=xk(:); yk=yk(:); fk=fk(:);
pk=exp(i*2*pi*xk);
qk=exp(i*2*pi*yk);
% basis function(interpoint distances)
basisfn=@(d) d;
% find coefficients
[ppk,pp]=meshgrid(pk,pk);
[qqk,qq]=meshgrid(qk,qk);
dist=sqrt(abs(pp-ppk).^2+abs(qq-qqk).^2);
a=[zeros(3,3) [ones(n,1) pk qk]'
      ones(n,1) pk qk basisfn(dist)];
abc=a\[zeros(3,1);fk];
conda=cond(a);
% evaluate at the following points
x=linspace(0,Lx,nx);
y=linspace(0,Ly,ny);
[y,x]=meshgrid(y,x);
p=exp(i*2*pi/Lx*x);
q=exp(i*2*pi/Ly*y);
[ppk,pp]=meshgrid(pk,p(:));
[qqk,qq]=meshgrid(qk,q(:));
dist=sqrt(abs(pp-ppk).^2+abs(qq-qqk).^2);
b=abc(1)+abc(2)*p(:)+abc(3)*q(:)+basisfn(dist)*abc(4:n+3);
b=reshape(real(b),nx,ny);

% define the indexes.
ii=1:nx;
jj=1:ny;
k=zeros(nx,ny);
k(:)=1:nx*ny;
kkk=[k k k;k k k;k k k];
kp0=kkk(nx+ii+1,ny+jj);
km0=kkk(nx+ii-1,ny+jj);
k0p=kkk(nx+ii,ny+jj+1);
k0m=kkk(nx+ii,ny+jj-1);
k=k(:);
kp0=kp0(:);
km0=km0(:);
k0p=k0p(:);
k0m=k0m(:);

% initial values
mp=0.4;
h=max(max(b))+mp-b; % the minimum shallow water depth is mp.
u=0*sqrt(0.992*sin(theta).*h/0.00293) ...
    +0.001*rand(nx,ny)+0*0.01*sin(2*pi/Lx*x);
v=zeros(nx,ny);
huv=[h(:)';u(:)';v(:)'];
huv=huv(:);

[ts huvs]=ode15s(@smagcurthree,[0 800],huv);

% surf the simulations
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
       axis([0 Lx 0 Ly 0 1.5]);
        title(['t=' num2str(ts(it))]);
    subplot(2,2,2),surf(x,y,b);
       xlabel('x');ylabel('y');zlabel('bed surface');
       %axis([0 Lx 0 Ly 0 1]);
    subplot(2,2,3),surf(x,y,u);
       xlabel('x');ylabel('y');zlabel('mean u');
        axis([0 Lx 0 Ly 0 1]);
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

return  
% plot the simulation
figure()
subplot(3,1.7,1);
contour(x,y,b);% contour the (meander) channel.
contour(x,y,u);% contour the mean lateral velocity.
contour(x,y,b-0.7+1,0.19,'k');% the shape of channel
contour(x,y,v);
contour(x,y,b-0.88+1,0.01,'k');
plot(y(21,12:29),h(21,12:29),'b');% plot depth, velocities at bends
plot(y(21,12:29),u(21,12:29),'b');
plot(y(21,12:29),v(21,12:29),'b');
qq=sqrt(huvs(:,2:3:end).^2+huvs(:,3:3:end).^2);% mean velocity
plot(ts,qq(:,200:500:end));%plot mean velocity history







% Simulations of turbulence flows by the Smagorinski model.
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
theta=0.01; % slope in x-direction
mag=0.9; % mid-depth of the channel
alpha=10; % centreline of the channel
beta=4; % width of the channel

% define the coordinate and indexes.
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

% define the channel
A=0; % straight channel zero A, and meandering channel nonzero A
f=A*cos(4*pi/Lx*x+pi/8);
b=-1+mag-mag*max(0,1-(y-f-alpha).^2/beta^2).^2;
 
% initial values
h=-b; % zero water level
u=0*sqrt(0.992*sin(theta).*h/0.00293) ...
    +0.01*rand(nx,ny)+0*0.01*sin(2*pi/Lx*x);
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
       axis([0 Lx 0 Ly -1 0.5]);
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







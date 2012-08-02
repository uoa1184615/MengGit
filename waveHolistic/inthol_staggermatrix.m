% Patches simulation of nonlinear problem. 
% MC and AJR, 18 Jul 2012

global n m r dx nu0 nu2 gam c1 theta g interOrd
tEnd=200; % the end time
nTime=200; % time grids
n=9;  % 1,5,9,... number of microscale grid points, not incl edges
m=8; % even number of patches
r=1/6 ;% ratio of patch size to macrogrid size
L=2*pi; % length of domain
D=L/m; % distance between neighbouring patches
d=2*r*D; % length of each patch
dx=d/(n+1); % patch size includes edges, to give microgrid size
nu0=0.01; % coefficient of bed drag
nu2=0.03; % coefficient of disspation
c1=1; % coefficient of nonlinear advection
theta=0.001; % the slope in the x-direction
g=1; % gravity
gam=1; % introduced coupling parameter
interOrd=3; % interpolation: 1=linear, 2=quadratic, and so on

% define the x-coords, omit patch edges
[hu,x]=meshgrid((0:m-1)*D,(0:n-1)*dx);
x=x+hu;

% initial values, hu is the correct size
k=4*pi/L;
hu(1:2:end)=0.2+0.1*sin(k*x(1:2:end));
hu(2:2:end)=   +0.1*sin(k*x(2:2:end));
hu=hu+0.02*(rand(size(hu))-0.5); % fun random perturbation

% do computation by defined function
ts=linspace(0,tEnd,nTime); 
hu0=hu(:);
[ts hus]=ode15s(@hol_staggermatrix,ts,hu);

% plot the simulation
for it=1:length(ts)
    figure(1)
    hu(:)=hus(it,:);
 plot(x(1:2:end,1:2:end),hu(1:2:end,1:2:end),'o' ...
        ,x(2:2:end,2:2:end),hu(2:2:end,2:2:end),'o' ...
        ,x(1:2:end,2:2:end),hu(1:2:end,2:2:end),'*' ...
        ,x(2:2:end,1:2:end),hu(2:2:end,1:2:end),'*')
      xlabel('x');ylabel('u and h');
    axis([0 L -0.2 0.4]);
    title(['t=' num2str(ts(it))]);
    drawnow
    %pause(0.1)
    M(:,it)=getframe(gcf);
end
%movie(M);
map=colormap; 
mpgwrite(M, map,'filename');


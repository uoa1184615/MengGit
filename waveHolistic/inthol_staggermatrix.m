global n m r dx nu0 nu2 gam c1 c2 theta g
n=9;  % 1,5,9,... number of microscale grid points, not incl edges
m=6; % even number of patches
r=1/6 ;% ratio of patch size to macrogrid size
L=2*pi; % length of domain
D=L/m;
d=2*r*D;
dx=d/(n+1); % patch size includes edges, to give microgrid size
nu0=0.1;
nu2=0.03;
c1=1;
c2=0.5;
theta=0.05;
g=1;
gam=1;

% define the x-coords, omit patch edges
[hu,x]=meshgrid((0:m-1)*D,(0:n-1)*dx);
x=x+hu;

% initial values, hu is the correct size
k=2*pi/L;
hu(1:2:end)=0.2+0.1*sin(k*x(1:2:end));
hu(2:2:end)=   +0.1*sin(k*x(2:2:end));

ts=linspace(0,160,320); 
hu0=hu(:);
[ts hus]=ode15s(@hol_staggermatrix,ts,hu);

% plot the simulation
for it=1:length(ts)
    figure(1)
    hu(:)=hus(it,:);
    plot(x(1:2:end,:),hu(1:2:end,:),'o' ...
        ,x(2:2:end,:),hu(2:2:end,:),'o')
      xlabel('x');ylabel('u and h');
    %axis([0 L+D -0.2 0.4]);
    title(['t=' num2str(ts(it))]);
    drawnow
    %pause(0.1)
    M(:,it)=getframe(gcf);
end
%movie(M);
map=colormap; 
mpgwrite(M, map,'filename');


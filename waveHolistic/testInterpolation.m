% Test patch interpolation formula
% AJR, 18 Jul 2012
n=5;  % 1,5,9,... number of microscale grid points, not incl edges
m=14; % even number of patches
r=1/6 ;% ratio of patch size to macrogrid size
L=2*pi; % length of domain
D=L/m; % distance between neighbouring patches
d=2*r*D; % length of each patch
dx=d/(n+1); % patch size includes edges, to give microgrid size
gam=1; % introduced coupling parameter
interOrd=3 % interpolation: 1=linear, 2=quadratic, and so on

% define the x-coords, do not omit patch edges
% as here I do want the patch edges
[X,x]=meshgrid((0:m-1)*D,(0:n+1)*dx);
x=x+X;

% initial values, hu is the correct size
hu=zeros(n,m);
% change a value to one, leaving all others zero
hu((n+1)/2,m/2)=1;

% now proceed with interpolation as in the function
% uu(i,j)=ith microgrid value in jth macropatch
uu=nan(n+2,m); 
i=2:n+1; j=1:m;
uu(i,j)=reshape(hu,n,m);

%periodic boundary condition from wrapping patch index
jp=[2:m, 1]; jm=[m, 1:m-1]; 
% sometimes it is more convenient to generate these as
% jpp=jp(jp(jp)); jppp=jp(jp(jpp)); and so on
jpp=[4:m, 1:3]; jmm=[m-2:m, 1:m-3]; 
jppp=[6:m,1:5]; jmmm=[m-4:m,1:m-5];

% coupling edge conditions on each tooth, h and u are the same
imid=(n+3)/2;
uubv=[gam*(uu(imid,jp)+uu(imid,jm))/2 
     +gam*r*(uu(imid,jp)-uu(imid,jm))/2 
    +(-1+r^2)*gam^2*(uu(imid,jpp)-uu(imid,jp)-uu(imid,jm)+uu(imid,jmm))/16
    +(-r+r^3)*gam^2*(uu(imid,jpp)-3*uu(imid,jp)+3*uu(imid,jm)-uu(imid,jmm))/48
    +(9-10*r^2+r^4)*gam^3*(uu(imid,jppp)+uu(imid,jpp)-2*uu(imid,jp)-2*uu(imid,jm)+uu(imid,jmm)+uu(imid,jmmm))/768
    +(9*r-10*r^3+r^5)*gam^3*(uu(imid,jppp)-5*uu(imid,jpp)+10*uu(imid,jp)-10*uu(imid,jm)+5*uu(imid,jmm)-uu(imid,jmmm))/3840
    ];
uu(1,j)=sum(uubv(1:2*interOrd,:));
uubv=[gam*(uu(imid,jp)+uu(imid,jm))/2 
     -gam*r*(uu(imid,jp)-uu(imid,jm))/2 
    +(-1+r^2)*gam^2*(uu(imid,jpp)-uu(imid,jp)-uu(imid,jm)+uu(imid,jmm))/16 
    -(-r+r^3)*gam^2*(uu(imid,jpp)-3*uu(imid,jp)+3*uu(imid,jm)-uu(imid,jmm))/48 
    +(9-10*r^2+r^4)*gam^3*(uu(imid,jppp)+uu(imid,jpp)-2*uu(imid,jp)-2*uu(imid,jm)+uu(imid,jmm)+uu(imid,jmmm))/768 
    -(9*r-10*r^3+r^5)*gam^3*(uu(imid,jppp)-5*uu(imid,jpp)+10*uu(imid,jp)-10*uu(imid,jm)+5*uu(imid,jmm)-uu(imid,jmmm))/3840
    ];
uu(n+2,j)=sum(uubv(1:2*interOrd,:))    

% plot the picture
    plot(x(1:2:end,:),uu(1:2:end,:),'o' ...
        ,x(2:2:end,:),uu(2:2:end,:),'*')
    xlabel('x');ylabel('u and h');

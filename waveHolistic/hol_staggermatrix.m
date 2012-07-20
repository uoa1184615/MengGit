% Patches simulation for nonlinear problem.
% MC and AJR, 18 Jul 2012

% define a function to contain the nonlinear problem.
function dhudt=hol_staggermatrix(t,hu) 
global n m r dx nu0 nu2 gam c1 theta g interOrd

% uu(i,j)=ith microgrid value in jth macropatch.
% uu(1,j) and uu(n+2,j) are the boundaries on every patch.
uu=nan(n+2,m); 
i=2:n+1; j=1:m;
uu(i,j)=reshape(hu,n,m);

% periodic boundary condition from wrapping patch index
jp=[2:m, 1]; jm=[m, 1:m-1]; 
% sometimes it is more convenient to generate these as
% jpp=jp(jp(jp)); jppp=jp(jp(jpp)); and so on
jpp=[4:m, 1:3]; jmm=[m-2:m, 1:m-3]; 
jppp=[6:m,1:5]; jmmm=[m-4:m,1:m-5];

% coupling edge conditions on each tooth, h and u are the same.
% linear interpolation when gam^2=0  
% quadratic interpolation when gam^3=0
imid=(n+3)/2;
uubv=[gam*(uu(imid,jp)+uu(imid,jm))/2 
     +gam*r*(uu(imid,jp)-uu(imid,jm))/2 
    +(-1+r^2)*gam^2*(uu(imid,jpp)-uu(imid,jp)-uu(imid,jm)+uu(imid,jmm))/16
    +(-r+r^3)*gam^2*(uu(imid,jpp)-3*uu(imid,jp)+3*uu(imid,jm)-uu(imid,jmm))/48
    +(9-10*r^2+r^4)*gam^3*(uu(imid,jppp)+uu(imid,jpp)-2*uu(imid,jp)-2*uu(imid,jm)+uu(imid,jmm)+uu(imid,jmmm))/768
    +(9*r-10*r^3+r^5)*gam^3*(uu(imid,jppp)-5*uu(imid,jpp)+10*uu(imid,jp)-10*uu(imid,jm)+5*uu(imid,jmm)-uu(imid,jmmm))/3840
    ];
uu(n+2,j)=sum(uubv(1:2*interOrd,:));
uubv=[gam*(uu(imid,jp)+uu(imid,jm))/2 
     -gam*r*(uu(imid,jp)-uu(imid,jm))/2 
    +(-1+r^2)*gam^2*(uu(imid,jpp)-uu(imid,jp)-uu(imid,jm)+uu(imid,jmm))/16 
    -(-r+r^3)*gam^2*(uu(imid,jpp)-3*uu(imid,jp)+3*uu(imid,jm)-uu(imid,jmm))/48 
    +(9-10*r^2+r^4)*gam^3*(uu(imid,jppp)+uu(imid,jpp)-2*uu(imid,jp)-2*uu(imid,jm)+uu(imid,jmm)+uu(imid,jmmm))/768 
    -(9*r-10*r^3+r^5)*gam^3*(uu(imid,jppp)-5*uu(imid,jpp)+10*uu(imid,jp)-10*uu(imid,jm)+5*uu(imid,jmm)-uu(imid,jmmm))/3840
    ];
uu(1,j)=sum(uubv(1:2*interOrd,:));    

% wave part of rhs is very easy by symmetry
dhudt(i,j)=-(uu(i+1,j)-uu(i-1,j))/(dx*2);

% other parts of equations require more craft: 
% h is when both same; u is when both different;
i0=2:2:n+1; i1=3:2:n+1; % even and odd within a patch
j0=1:2:m;   j1=2:2:m; % corresponding even and odd patches

% first some bed drag
dhudt(i1,j0)=dhudt(i1,j0)-nu0*uu(i1,j0);
dhudt(i0,j1)=dhudt(i0,j1)-nu0*uu(i0,j1);

% second some gravitational forcing
dhudt(i1,j0)=dhudt(i1,j0)+g*sin(theta);
dhudt(i0,j1)=dhudt(i0,j1)+g*sin(theta);

% now viscous drag, perhaps only works for n>=5 ??
dhudt(i1,j0)=dhudt(i1,j0) ...
    +nu2/(4*dx^2)*(uu(i1-2,j0)-2*uu(i1,j0)+uu(i1+2,j0));
uuxx=(uu(i0(1:end-2),j1)-2*uu(i0(2:end-1),j1)+uu(i0(3:end),j1))...
    /(4*dx^2);
dhudt(i0,j1)=dhudt(i0,j1) ...
    +nu2*[uuxx(1,:);uuxx;uuxx(end,:)];

% nonlinear advection terms
dhudt(i1,j0)=dhudt(i1,j0) ...
    -c1*uu(i1,j0).*(uu(i1+2,j0)-uu(i1-2,j0))/2/dx;
uux=(uu(i0(1:end-2),j1)-uu(i0(3:end),j1))/2/dx;
dhudt(i0,j1)=dhudt(i0,j1) ...
    -c1*uu(i0,j1).*[uux(1,:);uux;uux(end,:)];

% form time derivative into a vector
dhudt=reshape(dhudt(i,j),n*m,1);  

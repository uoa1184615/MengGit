function dhudt=hol_staggermatrix(t,hu) 
global n m r dx nu0 nu2 gam

% uu(i,j)=ith microgrid value in jth macropatch
%uu=nan(n+2,m); % linear interpolation
uu=nan(n+2,m); % quadratic interpolation
i=2:n+1; j=1:m;
uu(i,j)=reshape(hu,n,m);

%periodic boundary condition from wrapping patch index
jp=[2:m, 1]; jm=[m, 1:m-1]; % linear interpolation
jpp=[4:m, 1:3]; jmm=[m-2:m, 1:m-3]; % quadratic interpolation

% coupling edge conditions on each tooth, h and u are the same

%linear interpolation
imid=(n+3)/2;
%uu(1,j)  =((1+r)/2*uu(imid,jm)+(1-r)/2*uu(imid,jp))*gam ...
%                +(1-gam)*uu(imid,j);
%uu(n+2,j)=(1+r)/2*uu(imid,jp)+(1-r)/2*uu(imid,jm)*gam ...
%               +(1-gam)*uu(imid,j);

%quadratic interpolartion
uu(1,j)=gam*(uu(imid,jp)+uu(imid,jm))/2 ...
           +gam*r*(uu(imid,jp)-uu(imid,jm))/2 ...
           +(-1+r^2)*gam^2*(uu(imid,jpp)-uu(imid,jp) ...
              -uu(imid,jm)+uu(imid,jmm))/16 ...
            +(-r+r^3)*gam^2*(uu(imid,jpp)-3*uu(imid,jp) ...
              +3*uu(imid,jm)-uu(imid,jmm))/48;
uu(n+2,j)=gam*(uu(imid,jp)+uu(imid,jm))/2 ...
           -gam*r*(uu(imid,jp)-uu(imid,jm))/2 ...
           +(-1+r^2)*gam^2*(uu(imid,jpp)-uu(imid,jp) ...
              -uu(imid,jm)+uu(imid,jmm))/16 ...
            -(-r+r^3)*gam^2*(uu(imid,jpp)-3*uu(imid,jp) ...
              +3*uu(imid,jm)-uu(imid,jmm))/48;
        
% wave part of rhs is very easy by symmetry
dhudt(i,j)=-(uu(i+1,j)-uu(i-1,j))/(dx*2);

% other parts of equations require more craft: 
% h is when both same; u is when both different;
i0=2:2:n+1; i1=3:2:n+1; % even and odd within a patch
j0=1:2:m; j1=2:2:m; % corresponding even and odd patches

% first some bed drag
dhudt(i1,j0)=dhudt(i1,j0)-nu0*uu(i1,j0);
dhudt(i0,j1)=dhudt(i0,j1)-nu0*uu(i0,j1);

% now viscous drag, perhaps only works for n>=5 ??
dhudt(i1,j0)=dhudt(i1,j0)+nu2/(4*dx^2)*(uu(i1-2,j0) ...
                   -2*uu(i1,j0)+uu(i1+2,j0));
uuxx=uu(i0(1:end-2),j1)-2*uu(i0(2:end-1),j1)+uu(i0(3:end),j1);
dhudt(i0,j1)=dhudt(i0,j1)+nu2/(4*dx^2)*[uuxx(1,:);uuxx;uuxx(end,:)];

% form time derivative into a vector
dhudt=reshape(dhudt(i,j),n*m,1);  

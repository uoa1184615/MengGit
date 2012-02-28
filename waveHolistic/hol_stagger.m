function dhudt=hol_stagger(t,hu) 
global n m r dx nu0 nu2

% uu(i,j)=ith microgrid value in jth macropatch
uu=nan(n+2,m);
i=2:n+1; j=1:m;
uu(i,j)=reshape(hu,n,m);

%periodic boundary condition from wrapping patch index
jp=[2:m, 1]; jm=[m, 1:m-1];

% coupling edge conditions on each tooth, h and u are the same 
imid=(n+3)/2;
uu(1,j)  =(1+r)/2*uu(imid,jm)+(1-r)/2*uu(imid,jp);
uu(n+2,j)=(1+r)/2*uu(imid,jp)+(1-r)/2*uu(imid,jm);

% wave part of rhs is very easy by symmetry
dhudt(i,j)=-(uu(i+1,j)-uu(i-1,j))/(dx*2);

% other parts of equations require more craft: 
% h is when both same; u is when both different;
i0=2:2:n+1; i1=3:2:n+1; % even and odd within a patch
j0=1:2:m; j1=2:2:m; % corresponding even and odd patches

% first some bed drag
dhudt(i1,j0)-=nu0*uu(i1,j0);
dhudt(i0,j1)-=nu0*uu(i0,j1);

% now viscous drag, perhaps only works for n>=5 ??
dhudt(i1,j0)+=nu2/(4*dx^2)*(uu(i1-2,j0)-2*uu(i1,j0)+uu(i1+2,j0));
uuxx=uu(i0(1:end-2),j1)-2*uu(i0(2:end-1),j1)+uu(i0(3:end),j1);
dhudt(i0,j1)+=nu2/(4*dx^2)*[uuxx(1,:);uuxx;uuxx(end,:)];

% form time derivative into a vector
dhudt=reshape(dhudt(i,j),n*m,1);  

% multi-D linear spline prediction
xe=[0.1;0.2;0.5;0.3;0.45;0.6;0.8;0.62;0.9;1];
xk=zeros(3*size(xe,1),1);
xk(1:3:end)=xe;
xk(2:3:end)=xe;
xk(3:3:end)=xe;
ye1=[0.7;0.7;0.7;0.4;0.18;0.2;0.4;0.7;0.7;0.7];
ye2=[0.8;0.8;0.8;0.5;0.25;0.3;0.52;0.8;0.8;0.8];
ye3=[0.6;0.6;0.6;0.3;0.08;0.1;0.3;0.6;0.6;0.6];
yk=zeros(3*size(xe,1),1);
yk(1:3:end)=ye2;
yk(2:3:end)=ye1;
yk(3:3:end)=ye3;
fe1=[0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1];
fe2=[0.3;0.3;0.3;0.3;0.3;0.3;0.3;0.3;0.3;0.3];
fe3=[0.3;0.3;0.3;0.3;0.3;0.3;0.3;0.3;0.3;0.3];
fk=zeros(3*size(xe,1),1);
fk(1:3:end)=fe2;
fk(2:3:end)=fe1;
fk(3:3:end)=fe3;
% try this simpler version
xk=[.2 .2 .2 .2 .7 .7 .7 .7]
yk=[.1 .3 .6 .9 .1 .4 .6 .9]
fk=[0  1  0  0  0  0  1  0]
% force column vector and make periodic??
n=length(xk);
xk=xk(:); yk=yk(:); fk=fk(:);
pk=exp(i*2*pi*xk);
qk=exp(i*2*pi*yk);
% basis function(interpoint distances)
basisfn=@(d) log(d+1e-9).*d.^2;
% find coefficients
[ppk,pp]=meshgrid(pk,pk);
[qqk,qq]=meshgrid(qk,qk);
dist=sqrt(abs(pp-ppk).^2+abs(qq-qqk).^2);
a=[zeros(3,3) [ones(n,1) pk qk]'
      ones(n,1) pk qk basisfn(dist)];
abc=a\[zeros(3,1);fk];
conda=cond(a);
% evaluate at the following points
m=61;
[x,y]=meshgrid(linspace(0,1,m));
p=exp(i*2*pi*x);
q=exp(i*2*pi*y);
[ppk,pp]=meshgrid(pk,p(:));
[qqk,qq]=meshgrid(qk,q(:));
dist=sqrt(abs(pp-ppk).^2+abs(qq-qqk).^2);
z=abc(1)+abc(2)*p(:)+abc(3)*q(:)+basisfn(dist)*abc(4:n+3);
z=reshape(real(z),m,m);
clf()
%subplot(1.7,1.7,1);
%surf(x,y,z);
contour(x,y,z), hold on
xlabel('x'), ylabel('y'), zlabel('f')
plot(xk,yk,'+'), hold off
%plot3(x(1,:),y(:,1),real(z)','theta=65,alpha=20');
%cf=gcf();
%cf.color_map=rgbmap(32);

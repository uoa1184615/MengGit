%% Construct a river-like channel
%% MC and AJR, 26 JUly 2012
% multi-D linear spline prediction
xk=[.2 .2 .2 .2 .5 .5 .5 .5 .3 .3 .3 .3 .4 .4 .4 .4 .7 .7 .7 .7 .8 .8 .8 .8 .9 .9 .9 .9];
yk=[.1 .3 .6 .9 .3 .5 .7 .9 .1 .4 .6 .9 .1 .4 .6 .9 .1 .4 .6 .9 .1 .3 .6 .9 .1 .3 .6 .9];
fk=[1  0  1  1  1  1  0  1 1 0 1 1 1 1 0 1 1 1 0 1 1 1 0 1 1 0 1 1];
% force column vector and make periodic??
n=length(xk);
xk=xk(:); yk=yk(:); fk=fk(:);
pk=exp(i*2*pi*xk);
qk=exp(i*2*pi*yk);
% find coefficients
[ppk,pp]=meshgrid(pk,pk);
[qqk,qq]=meshgrid(qk,qk);
dist=sqrt(abs(pp-ppk).^2+abs(qq-qqk).^2);
a=[zeros(3,3) [ones(n,1) pk qk]'
      ones(n,1) pk qk dist.^3];
abc=a\[zeros(3,1);fk];
conda=cond(a);
% evaluate at the following points
%m=61;
x=linspace(0,40,40);
y=linspace(0,20,40);
[y,x]=meshgrid(y,x);
p=exp(i*2*pi/40*x);
q=exp(i*2*pi/20*y);
[ppk,pp]=meshgrid(pk,p(:));
[qqk,qq]=meshgrid(qk,q(:));
dist=sqrt(abs(pp-ppk).^2+abs(qq-qqk).^2);
z=abc(1)+abc(2)*p(:)+abc(3)*q(:)+(dist.^3)*abc(4:n+3);
z=reshape(real(z),40,40);
clf()
%subplot(1.7,1.7,1);
%surf(x,y,z);
contour(x,y,z), hold on
xlabel('x'), ylabel('y'), zlabel('f')
plot(xk,yk,'+'), hold off
%plot3(x(1,:),y(:,1),real(z)','theta=65,alpha=20');
%cf=gcf();
%cf.color_map=rgbmap(32);

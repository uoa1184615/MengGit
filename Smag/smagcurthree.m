% Simulations of turbulence flows by the Smagorinski model.
% MC and AJR, 19 July 2012
% define a function include the governing equations.
function dhuv=smagcurthree(t,huv)
global k kp0 km0 k0p k0m dx dy theta g b gam

% create a new vector
h=huv(1:3:end);
u=huv(2:3:end);
v=huv(3:3:end);

dhdt=-h(k).*(u(kp0)-u(km0))/dx/2-u(k).*(h(kp0)-h(km0))/dx/2 ...
         -h(k).*(v(k0p)-v(k0m))/dy/2-v(k).*(h(k0p)-h(k0m))/dy/2;
  
 qq=sqrt(u(k).^2+v(k).^2);
 rqq=1./sqrt(u(k).^2+v(k).^2);
     
dudt=-0.00293*gam*u(k).*qq./h(k) ...
     +0.992*g*sin(theta) ...
     -0.992*g*cos(theta)*(h(kp0)-h(km0))/dx/2 ...
     -0.992*g*cos(theta)*(b(kp0)-b(km0))/dx/2 ...
     -1.04*v(k).*(u(k0p)-u(k0m))/dy/2 ...
     -1.025*u(k).*(u(kp0)-u(km0))/dx/2 ...
     -0.0115*u(k).*(v(k0p)-v(k0m))/dy/2 ...
     +0.237*h(k).*qq.*(u(k0p)-2*u(k)+u(k0m))/dy^2 ...
     +0.0266*h(k).*qq.*(u(kp0)-2*u(k)+u(km0))/dx^2 ...
     -0.00784*h(k).*u(k).*v(k).*rqq.*(v(k0p)-2*v(k)+v(k0m))/dy^2 ...
     +0.00119*h(k).*u(k).*v(k).*rqq.*(v(kp0)-2*v(k)+v(km0))/dx^2 ...
     +0.468*h(k).*u(k).^2.*rqq.*(u(kp0)-2*u(k)+u(km0))/dx^2;
   
     
 dvdt=-0.00293*gam*v(k).*qq./h(k) ...
     -0.992*g*cos(theta)*(h(k0p)-h(k0m))/dy/2 ...
      -0.992*g*cos(theta)*(b(k0p)-b(k0m))/dy/2 ...
      -1.04*v(k).*(v(k0p)-v(k0m))/dy/2 ...
       -1.025*u(k).*(v(kp0)-v(km0))/dx/2 ...
       -0.0084*v(k).*(u(kp0)-u(km0))/dx/2 ...
       +0.2488*h(k).*qq.*(v(k0p)-2*v(k)+v(k0m))/dy^2 ...
       +0.00606*h(k).*qq.*(v(kp0)-2*v(k)+v(km0))/dx^2 ...
       -0.0119*h(k).*u(k).^2.*rqq.*(v(kp0)-2*v(k)+v(km0))/dx^2 ...
       -0.243*h(k).*u(k).*v(k).*rqq.*(u(k0p)-2*u(k)+u(k0m))/dy^2 ...
       +0.235*h(k).*u(k).*v(k).*rqq.*(u(kp0)-2*u(k)+u(km0))/dx^2;

         
dhuv=[dhdt';dudt';dvdt'];
dhuv=dhuv(:);

function dhuv=smagcurthree(t,huv)
global k kp0 km0 k0p k0m dx dy theta g b c
h=huv(1:3:end);
u=huv(2:3:end);
v=huv(3:3:end);

dhdt=-h(k).*(u(kp0)-u(km0))/dx/2-u(k).*(h(kp0)-h(km0))/dx/2 ...
         -h(k).*(v(k0p)-v(k0m))/dy/2-v(k).*(h(k0p)-h(k0m))/dy/2;
     
 dudt=-0.0029*u(k).*sqrt(u(k).^2+v(k).^2)./h(k) ...
         -1.03*v(k).*(u(k0p)-u(k0m))/dy/2 ...
         -1.045*u(k).*(u(kp0)-u(km0))/dx/2+0.985*g*sin(theta) ...
         -0.985*g*cos(theta)*(h(kp0)-h(km0))/dx/2 ...
         -0.985*g*cos(theta)*(b(kp0)-b(km0))/dx/2 ...
         +c*(u(kp0)-2*u(k)+u(km0))/dx^2 ...
         +c*(u(k0p)-2*u(k)+u(k0m))/dy^2;    
     
 dvdt=-0.0029*v(k).*sqrt(u(k).^2+v(k).^2)./h(k) ...
         -1.03*u(k).*(v(kp0)-v(km0))/dx/2 ...
         -1.045*v(k).*(v(k0p)-v(k0m))/dy/2 ...
         -0.989*g*cos(theta)*(h(k0p)-h(k0m))/dy/2 ...
         -0.989*g*cos(theta)*(b(k0p)-b(k0m))/dy/2 ...
         +c*(v(kp0)-2*v(k)+v(km0))/dx^2 ...
         +c*(v(k0p)-2*v(k)+v(k0m))/dy^2;
         
dhuv=[dhdt';dudt';dvdt'];
dhuv=dhuv(:);
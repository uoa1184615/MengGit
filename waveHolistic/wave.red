write "try to develop good numerics of wave-like PDEs using
a staggered element approach.  It is a bit hairy.  MC & AJR
Oct 2011";
on div; off allfac; on revpri;
let df(sign(~x),~y)=>0;
write "xi=(x-X_j)/dx is internal variable, |xi|<1/2";
depend xi,x;  let df(xi,x)=>1/dx;
write "amplitudes are defined at ??";
operator hh; operator uup;
depend hh,t; depend uup,t;
let { df(hh(~k),t)=>sub(j=k,gh)
    , df(uup(~k),t)=>sub(j=k,gup)
    };
write "linear approximation";
h:=hh(j); up:=uup(j);
gh:=gup:=0;
let gam^2=>0;
write "compute some residuals, but shifted so that the xi
variables are the same for the terms";
reshr:=sub(xi=xi+1/2,df(h,t))+df(up,x);
resur:=sub(xi=xi-1/2,df(up,t))+df(h,x);
reshl:=sub({xi=xi-1/2,j=j+1},df(h,t))+df(up,x);
resul:=sub({xi=xi+1/2,j=j-1},df(up,t))+df(h,x);
reshc:=sub(xi=1/2,h)-sub(xi=-1/2,h)-gam/2*(
     sub(xi=+1/2,h)+sub({xi=-1/2,j=j+1},h)
    -sub(xi=-1/2,h)-sub({xi=+1/2,j=j-1},h));
resuc:=sub(xi=1/2,up)-sub(xi=-1/2,up)-gam/2*(
     sub(xi=+1/2,up)+sub({xi=-1/2,j=j+1},up)
    -sub(xi=-1/2,up)-sub({xi=+1/2,j=j-1},up));
end;
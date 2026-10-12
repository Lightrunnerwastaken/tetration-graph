/* Real-only finite Koenigs approximants, no complex continuation asserted. */
hiter(r,z,k)={for(j=1,k,z=hval(r,z));return(z);}
ha=0;hla=0;hlaerr=0;
haready(r)={
  my(y=1,y2,hh=1e-5,d1,d2);
  hinit(r-1);
  for(j=1,160,y2=hval(r-1,y);if(abs(y2-y)<1e-25,y=y2;break);y=y2);
  ha=y;
  d1=(hval(r-1,y+hh)-hval(r-1,y-hh))/(2*hh);
  d2=(hval(r-1,y+hh/2)-hval(r-1,y-hh/2))/hh;
  hla=(4*d2-d1)/3;hlaerr=abs(d2-d1);
  if(!(0<hla && hla<1),error("haready: not attracting"));
  return(ha);
}
hfinite(r,t,k)={
  my(a=1,b=sexp(1),m,goal=ha+hla^t*(hiter(r-1,1,k)-ha),v);
  if(t==0,return(1));if(t==1,return(sexp(1)));
  if(!(0<t && t<1),error("hfinite: height outside [0,1]"));
  for(j=1,42,m=(a+b)/2;v=hiter(r-1,m,k);if(v<goal,a=m,b=m));
  return((a+b)/2);
}
hreal(r,x,k)={
  my(n=floor(x),t=x-n,y);
  if(n<0 || n>2,error("hreal: numerical demonstration only on [0,2]"));
  y=hfinite(r,t,k);return(hiter(r-1,y,n));
}
hback(r,z,k)={for(j=1,k,z=hinv(r,z));return(z);}
hrealback(r,t,k)={
  my(z=ha+hla^t*(hiter(r-1,1,k)-ha));
  return(hback(r-1,z,k));
}
halflast=0;
hrealset(r,t,k)={halflast=hrealback(r,t,k);return(0);}

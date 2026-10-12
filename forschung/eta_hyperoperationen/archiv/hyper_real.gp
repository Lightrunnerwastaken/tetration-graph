/* Own uniform hyperoperation interface. Upstream engine is a local snapshot.
   Odd ranks >=7: regular iteration at left real repelling fixed point.
   Even ranks >=8: regular iteration at a positive attracting fixed point.
   Rank 8 additionally has an experimental complex two-arc Abel selection.
   All estimates are numerical, not interval certificates. */
hrp=vector(16); hrlam=vector(16); hrpsi=vector(16); hrsig=vector(16);
hrrad=vector(16); hrc=vector(16); hrready=vector(16); hrerror=vector(16);
hrtarget=vector(16); hrdepth=vector(16); hrshift=vector(16);
hrkind=vector(16);
hrmaxdigits=vector(16);hrdebug=0;

hval(r,z,extra=0)={
  my(w,y,k=0,la);
  if(r==3,return(sexp(1)^z));
  if(r==4,return(sexp(z)));
  if(r==5,return(pent(z,extra)));
  if(r==6,return(hex(z)));
  if(!hrready[r],error("hval: rank not initialized"));
  if(hrkind[r]==2,return(o8knhex(z)));
  la=hrlam[r];w=hrc[r]*la^z;
  if(la>1,
    while(abs(w)>hrrad[r],w/=la;k++;if(k>400,error("hval: depth >400"))),
    while(abs(w)>hrrad[r],w*=la;k++;if(k>400,error("hval: depth >400"))));
  for(j=1,extra,if(la>1,w/=la,w*=la);k++);
  y=hrp[r]+subst(hrpsi[r],x,w);
  if(la>1,for(j=1,k,y=hval(r-1,y)),for(j=1,k,y=hinv(r-1,y)));
  return(levelchk(y,"hval"));
}
hinv(r,z)={
  my(y=z,k=0,la,w);
  if(r==3,return(log(z)/log(sexp(1))));
  if(r==4,return(slog(z)));
  if(r==5,return(ipent(z)));
  if(r==6,return(ihex(z)));
  if(!hrready[r],error("hinv: rank not initialized"));
  if(hrkind[r]==2,return(o8knihex(z)));
  la=hrlam[r];
  while(abs(y-hrp[r])>hrrad[r],
    if(la>1,y=hinv(r-1,y),y=hval(r-1,y));
    k++;if(k>400,error("hinv: depth >400")));
  w=subst(hrsig[r],x,y-hrp[r]);
  return(log(w/hrc[r])/log(la)+if(la>1,k,-k));
}
hinit(r,sr=1/8,rr=1/100)={
  my(y=1,y2,sc,res,ps,eps,la,err=0,target,k=0,w,c0,R,nmax,settled=0);
  if(r<=4,return(precis));
  if(r==5,return(pentensure()));
  if(r==6,return(hexensure()));
  if(hrready[r],return(hrerror[r]));
  if(hrdebug,print("hinit: rank ",r," begin"));
  hinit(r-1);
  target=max(8,min(precis-10,if(r==7,floor(hexdigits)-6,floor(hrerror[r-1])-5)));
  if(hrmaxdigits[r]>0,target=min(target,hrmaxdigits[r]));
  eps=10.0^(-target-3);hrtarget[r]=target;
  for(j=1,300,
    /* Numerical branch selection in the investigated range: the measured
       attracting real fixed points lie below 4. Stop the exploratory orbit
       before enormous higher iterates. This is not a global absence proof. */
    if(!(r%2) && y>=4,return(hcomplexinit(r)));
    y2=if(r%2,hinv(r-1,y),iferr(hval(r-1,y),E,10^100));
    if(abs(imag(y2))>1e-10,
      if(!(r%2),return(hcomplexinit(r)),error("hinit: no real fixed point on selected branch")));
    y2=real(y2);
    if(abs(y2)>30,
      if(!(r%2),return(hcomplexinit(r)),error("hinit: no real fixed point on selected branch")));
    if(abs(y2-y)<1e-9,y=y2;settled=1;break);y=y2);
  if(!settled,error("hinit: fixed-point orbit failed"));
  if(hrdebug,print("hinit: rank ",r," candidate ",precision(y,10)," target ",target));
  if(r%2,res=pentfixser(t->hinv(r-1,t),y,sr,1,eps,"hinit odd"),
         res=pentfixser(t->hval(r-1,t),y,sr,1,eps,"hinit even"));
  hrp[r]=res[1];sc=res[2];R=res[3];
  if(hrdebug,print("hinit: rank ",r," local series ready"));
  la=polcoeff(sc,1);
  if(r%2,la=1/la);
  if((r%2 && la<=1)||(!(r%2) && (la<=0 || la>=1)),error("hinit: wrong multiplier"));
  hrlam[r]=la;
  rr=min(rr,R/20);
  ps=pentpsi(sc,rr,eps,poldegree(sc)-1);
  if(!ps[2],error("hinit: local Poincare series failed"));
  hrpsi[r]=ps[1];hrsig[r]=Pol(serreverse(ps[1]+O(x^(poldegree(ps[1])+1))));
  hrrad[r]=rr;
  y=1;
  while(abs(y-hrp[r])>rr,
    y=if(la>1,hinv(r-1,y),hval(r-1,y));k++;
    if(k>400,error("hinit: normalization orbit >400")));
  w=subst(hrsig[r],x,y-hrp[r]);
  c0=w/la^if(la>1,-k,k);
  hrc[r]=c0;hrdepth[r]=k;hrready[r]=1;hrkind[r]=1;
  for(j=0,3,
    y=j/4;
    err=max(err,abs(hval(r,y)-hval(r,y,2))/(1+abs(hval(r,y))));
    err=max(err,abs(hinv(r,hval(r,y))-y)));
  err=max(err,max(abs(hval(r,0)-1),abs(hval(r,1)-sexp(1))));
  hrerror[r]=if(err==0,target,min(target,-log(err)/log(10)));
  return(hrerror[r]);
}
hprobe()={my(y=1);for(j=1,20,y=ihex(y));return(y);}
hminimum(r,z,h)={
  my(a=hval(r,z-h),b=hval(r,z),c=hval(r,z+h));
  return([b-z,(c-a)/(2*h)-1,(c-2*b+a)/h^2]);
}
hjet(r,z,rad,N)={
  my(P=pentcauchy(t->hval(r,t),z,rad,N));
  return([polcoeff(P,0)-z,polcoeff(P,1)-1,2*polcoeff(P,2)]);
}
hlastjet=0;
hjetset(r,z,rad,N)={hlastjet=hjet(r,z,rad,N);return(0);}

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

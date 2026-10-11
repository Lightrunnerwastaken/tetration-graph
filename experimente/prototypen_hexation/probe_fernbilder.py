"""Versuch zur e-Grenze: knsolve ohne Funktionalgleichungen, deren Bildpunkt weiter als gmax·Bogenradius vom Mittelpunkt liegt.
 kn2solve ist eine zur Laufzeit definierte Kopie von knsolve (Fork-Datei unverändert); Fehler-Vektor wie in knhexinit (zwei Bögen, 7 Testpunkte).
 Aufruf: probe_fernbilder.py basis dps knforce gmax1,gmax2,...
"""
import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe

basis, dps, kraft, gmaxs = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4].split(",")
mp.mp.dps = dps + 10
t0 = time.time()
KN2 = ("if(1, kn2solve=((rad,K,N,hcut,hcut2,gmax)->if(1, my(pts,W,eqs,U,rows,scl,Mt,bv,cn,sol,A,cc,rho,tp,u1,u2,uk1,uk2,val,e,nx,nk,g,res,nth=0,nfe=0,nweg=0,sav); "
       "pts=knarc(rad,ceil(0.8*(K+1+2*(N+1)))); W=knwalk(pts); eqs=List(); rho=rad; "
       "for(j=1,#pts, if(kns*imag(W[j])>=hcut, listput(eqs,[1,pts[j],W[j]]); nth++); "
       "if(kns*imag(W[j])<hcut2, if(real(pts[j])>knc, g=ipent(pts[j]); tp=1, g=pent(pts[j]); tp=-1); "
       "if(abs(g-knc)<=gmax*rad, listput(eqs,[2,pts[j],g,tp]); rho=max(rho,abs(g-knc)); nfe++, nweg++))); "
       "U=K+1+2*(N+1); rows=2*#eqs+1; "
       "scl=vector(N+1,n,vecmax(concat([1e-300],vector(#eqs,j,if(eqs[j][1]==1,abs(knen(eqs[j][3],n-1)),0))))); "
       "Mt=matrix(rows,U); bv=vector(rows)~; "
       "for(j=1,#eqs, tp=eqs[j][1]; u1=(eqs[j][2]-knc)/rho; uk1=1; "
       "if(tp==1, for(k=1,K+1, Mt[2*j-1,k]=real(uk1); Mt[2*j,k]=imag(uk1); uk1*=u1); "
       "for(n=1,N+1, e=knen(eqs[j][3],n-1)/scl[n]; Mt[2*j-1,K+1+2*n-1]=-real(e); Mt[2*j,K+1+2*n-1]=-imag(e); "
       "Mt[2*j-1,K+1+2*n]=imag(e); Mt[2*j,K+1+2*n]=-real(e)); bv[2*j-1]=real(eqs[j][3]); bv[2*j]=imag(eqs[j][3]), "
       "u2=(eqs[j][3]-knc)/rho; uk2=1; for(k=1,K+1, val=uk1-uk2; Mt[2*j-1,k]=real(val); Mt[2*j,k]=imag(val); uk1*=u1; uk2*=u2); "
       "bv[2*j-1]=eqs[j][4]; bv[2*j]=0)); "
       "nx=1; nk=0; while((abs(nx-knc)>0.9*rad)&&(nk<6), nx=ipent(nx); nk++); "
       "u1=(real(nx)-knc)/rho; for(k=1,K+1, Mt[rows,k]=u1^(k-1)); bv[rows]=-nk; "
       "cn=vector(U,k,vecmax(abs(Mt[,k]))+1e-300); Mt=matrix(rows,U,i,k,Mt[i,k]/cn[k]); Mt[rows,]=Mt[rows,]*1e6; bv[rows]*=1e6; "
       "sav=default(realprecision); default(realprecision,2*precis); sol=knlsq(Mt,bv); default(realprecision,sav); "
       "sol=vector(U,k,precision(sol[k]/cn[k],precis)); A=sum(k=1,K+1,sol[k]*x^(k-1)); cc=vector(N+1,n,(sol[K+1+2*n-1]+I*sol[K+1+2*n])/scl[n]); "
       "res=0; for(j=1,#eqs, tp=eqs[j][1]; val=if(tp==1, subst(A,x,(eqs[j][2]-knc)/rho)-eqs[j][3]-sum(n=1,N+1,cc[n]*knen(eqs[j][3],n-1)), "
       "subst(A,x,(eqs[j][2]-knc)/rho)-subst(A,x,(eqs[j][3]-knc)/rho)-eqs[j][4]); res=max(res,abs(val))); "
       "[A,cc,rho,res,nth,nfe,nweg])); 0)")
VERS = ("if(1, kr1=knpar[3]*knr; kp1=knplan(kr1,precis-3); kr2=kr1+0.05*knr; kp2=knplan(kr2,precis-3); "
        "kT1=kn2solve(kr1,knpar[1],knpar[2],kp1[5]+0.6*(kp1[6]-kp1[5]),kp1[5]+0.85*(kp1[6]-kp1[5]),{g}); "
        "kT2=kn2solve(kr2,knpar[1],knpar[2],kp2[5]+0.6*(kp2[6]-kp2[5]),kp2[5]+0.85*(kp2[6]-kp2[5]),{g}); "
        "ktst=concat(vector(5,j,knc+(j-3)/2*0.9*min(kr1,kr2)),[knc+0.5*kr1*exp(Pi*I/3),knc+0.5*kr1*exp(2*Pi*I/3)]); "
        "kdd=0; for(j=1,#ktst, kdd=max(kdd,abs(subst(kT1[1],x,(ktst[j]-knc)/kT1[3])-subst(kT2[1],x,(ktst[j]-knc)/kT2[3])))); "
        "-log(max(kdd,kT1[4]))/log(10))")
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=dps)
try:
    k = gp.eval_batch(basis, [f"if(1, knforce={kraft}; knmaxk=2000; 0)", "hexinit(sexp(1),0,1)", KN2, "knpar[1]", "knpar[2]", "knc"])
    print(f"Basis {basis}, dps {dps}, knforce {kraft}: Original {mp.nstr(k[1], 5)} Stellen, K = {int(k[3].real)}, N = {int(k[4].real)}, "
          f"knc = {mp.nstr(k[5], 8)} [{time.time() - t0:.0f} s]", flush=True)
    for g in gmaxs:
        w = gp.eval_batch(basis, [VERS.format(g=g), "kT1[4]", "kT2[4]", "kT1[6]", "kT1[7]", "kT2[7]", "kT1[3]/kr1"])
        print(f"  gmax = {g:5s}: {mp.nstr(w[0], 5)} Stellen; Reste {mp.nstr(w[1], 3)} / {mp.nstr(w[2], 3)}; FE-Gleichungen {int(w[3].real)}, "
              f"weggelassen {int(w[4].real)} / {int(w[5].real)}; ρ/Bogen {mp.nstr(w[6], 4)} [{time.time() - t0:.0f} s]", flush=True)
finally:
    gp.close()

"""E9: Kneser-Grössen für grosse Basen: T_b'(0), T_b''(0), Wendepunkt x_w(b) (Taylor um 0 und um −0.5, Newton auf T'').
Aufruf: python e9_grosse_basis.py DPS b1 b2 ...  (b als GP-Ausdruck)"""
import sys, json, mpmath as mp
from motor import werte, schliessen, HIER
dps=int(sys.argv[1]); mp.mp.dps=dps
K=32; R=mp.mpf("0.2"); out={}
def taylor(b,z0):
    zs=[z0+R*mp.expjpi(2*mp.mpf(j)/K) for j in range(K)]
    v=werte(b,[f"sexp({mp.nstr(mp.re(z),dps+5)}+({mp.nstr(mp.im(z),dps+5)})*I)" for z in zs],dps)
    return [mp.re(mp.fsum(v[j]*mp.expjpi(-2*mp.mpf(j*k)/K) for j in range(K))/K/R**k) for k in range(12)]
for b in sys.argv[2:]:
    try:
        a=taylor(b,mp.mpf(0))
        # Wendepunkt: Nullstelle von T'' = Σ k(k-1) a_k h^{k-2}, Newton ab h=0, ggf. neu zentrieren
        c=mp.mpf(0); A=a
        for it in range(6):
            p2=lambda h: mp.fsum(k*(k-1)*A[k]*h**(k-2) for k in range(2,12))
            h=mp.findroot(p2, mp.mpf(0))
            if abs(h)<0.05: break
            c=c+h; A=taylor(b,c)
        xw=c+h
        lb=mp.log(mp.mpf(eval(b.replace('^','**'))) if 'exp' not in b else mp.mpf(1))
        print(f"b={b}: T'(0)={mp.nstr(a[1],15)} T''(0)={mp.nstr(2*a[2],12)} x_w={mp.nstr(xw,15)}", flush=True)
        out[b]=[mp.nstr(a[1],20),mp.nstr(2*a[2],20),mp.nstr(xw,20)]
    except Exception as e:
        print(b,"FEHLER",str(e)[-150:], flush=True)
    schliessen()
(HIER/"daten").mkdir(exist_ok=True)
(HIER/"daten"/f"e9_{dps}.json").write_text(json.dumps(out,indent=1))

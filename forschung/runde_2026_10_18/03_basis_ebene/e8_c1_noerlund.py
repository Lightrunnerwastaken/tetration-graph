"""E8: Erste Basisableitung c1(x) = ∂_b T_b(x) bei η als Nörlund-Hauptlösung.
c1(x) = T_η'(x) G(x),  G(x+1) − G(x) = r(x) := (e/η) T_η(x)/T_η'(x),  G(0) = 0,  G = Hauptlösung (ohne periodischen Anteil).
G(x) = lim_N [ Σ_{n<N} (r(n) − r(x+n)) + ∫_N^{N+x} r − (r(N+x)−r(N))/2 + (r'(N+x)−r'(N))/12 − (r'''(N+x)−r'''(N))/720 ].
Aufruf: python e8_c1_noerlund.py N x1 x2 ..."""
import sys, mpmath as mp
from motor import werte, schliessen
mp.mp.dps=40; DP=40
B="(exp(exp(-1))+1/100)"   # Motor-Sitzung nur für sexpeta (hängt nicht von der Basis ab)
eta=mp.exp(mp.exp(-1)); E=mp.e
K=24; RR=mp.mpf("0.25")
def taylor(z0):
    zs=[z0+RR*mp.expjpi(2*mp.mpf(j)/K) for j in range(K)]
    v=werte(B,[f"sexpeta({mp.nstr(mp.re(z),45)}+({mp.nstr(mp.im(z),45)})*I)" for z in zs],DP)
    return [mp.fsum(v[j]*mp.expjpi(-2*mp.mpf(j*k)/K) for j in range(K))/K/RR**k for k in range(8)]
def rjet(t):
    """r(t), r'(t), r'''(t) aus Taylor von T um t."""
    a=taylor(mp.mpf(t)); a=[mp.re(z) for z in a]
    # r = (e/eta) T/T'  als Reihe in h
    T=a[:6]; Tp=[(k+1)*a[k+1] for k in range(6)]
    # Division T/Tp als Potenzreihe
    q=[mp.mpf(0)]*5
    for n in range(5):
        q[n]=(T[n]-mp.fsum(q[k]*Tp[n-k] for k in range(n)))/Tp[0]
    f=E/eta
    return f*q[0], f*q[1], f*6*q[3]
def rint(t0,x):
    # ∫_{t0}^{t0+x} r : Gauss-Legendre mit 8 Knoten über Taylor um Mitte
    m=t0+x/2; a=[mp.re(z) for z in taylor(mp.mpf(m))]
    T=a[:8]; Tp=[(k+1)*a[k+1] for k in range(7)]+[0]
    q=[mp.mpf(0)]*7
    for n in range(7):
        q[n]=(T[n]-mp.fsum(q[k]*Tp[n-k] for k in range(n)))/Tp[0]
    h=mp.mpf(x)/2
    return (E/eta)*mp.fsum(q[k]*(h**(k+1)-(-h)**(k+1))/(k+1) for k in range(7))
N=int(sys.argv[1])
for xs in sys.argv[2:]:
    x=mp.mpf(xs)
    S=mp.fsum(rjet(n)[0]-rjet(n+x)[0] for n in range(N))
    r1,d1,t1=rjet(N+x); r0,d0,t0=rjet(N)
    G=S+rint(N,x)-(r1-r0)/2+(d1-d0)/12-(t1-t0)/720
    Tp=mp.re(taylor(x)[1])
    print(f"x={xs} N={N}: G={mp.nstr(G,15)}  c1=T'G={mp.nstr(Tp*G,15)}", flush=True)
schliessen()

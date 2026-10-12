"""k14 — Test zu R2 (Renormierung der Kopfprofile): lineare Abbildung Δ_n ↦ Δ_{2n} (Δ_n = S_{2n} − S_n) aus Gs Profilen (spiral_analysis.json, nur gelesen;
vorher als profiles_kub.npz im Notizordner extrahiert: n, S (40 Koeffizienten), B). Eigenwerte vs. 2^{−a}, 2^{−1}, 2^{−2a}."""
import numpy as np
D=np.load(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'k14_profile_G_kubisch.npz')); n=list(D['n']); S=D['S']; B=D['B']
L=0.31813150520476413+1.3372357014306895j; mu=np.sqrt(L)
a=1/(2j*np.pi/(2j*np.pi-np.log(mu)))
print('erwartet 2^-a = %s (|.|=%.4f), 2^-1 = 0.5, 2^-2a = %s'%(np.round(2**-a,4),abs(2**-a),np.round(2**(-2*a),4)))
K=16
def s(m): return S[n.index(m)][1:K]
def Dl(m): return s(2*m)-s(m)
for pairs,name in [([(50,100),(100,200),(200,400),(150,300)],'alle 4 Paare'),([(100,200),(200,400),(150,300)],'n>=100: 3 Paare'),([(100,200),(200,400)],'Kette 100-200-400'),([(150,300),(200,400)],'n>=150')]:
    D0=np.array([Dl(p) for p,q in pairs]).T; D1=np.array([Dl(q) for p,q in pairs]).T
    T=np.linalg.lstsq(D0,D1,rcond=None)[0]
    ev=np.linalg.eigvals(T)
    res=np.linalg.norm(D0@T-D1)/np.linalg.norm(D1)
    print('%-18s Eigenwerte: %s   (|.|: %s)  rel.Rest %.1e'%(name, ', '.join('%.4f%+.4fi'%(e.real,e.imag) for e in ev), ', '.join('%.4f'%abs(e) for e in ev), res))
# Skalarer Test: Verhältnis aufeinanderfolgender Differenzen in jeder Komponente
for k in [1,2,3]:
    r1=(s(400)[k-1]-s(200)[k-1])/(s(200)[k-1]-s(100)[k-1]); r2=(s(800)[k-1]-s(400)[k-1])/(s(400)[k-1]-s(200)[k-1]); r3=(s(600)[k-1]-s(300)[k-1])/(s(300)[k-1]-s(150)[k-1])
    print('Koeff %d: Differenzenquotient 100-200-400: %.4f%+.4fi, 200-400-800: %.4f%+.4fi, 150-300-600: %.4f%+.4fi'%(k,r1.real,r1.imag,r2.real,r2.imag,r3.real,r3.imag))
# Amplitude
def b(m): return B[n.index(m)]
print('B: Differenzenquotienten 100-200-400: %.4f, 200-400-800: %.4f, 150-300-600: %.4f'%((b(400)-b(200))/(b(200)-b(100)),(b(800)-b(400))/(b(400)-b(200)),(b(600)-b(300))/(b(300)-b(150))))

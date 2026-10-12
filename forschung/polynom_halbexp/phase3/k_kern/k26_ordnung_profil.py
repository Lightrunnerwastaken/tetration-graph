"""k26 — Wachstumsordnung von Gs Eckprofil S_1000 (spiral_analysis.json, nur gelesen; Auszug S1000.npy aus dem Notizordner):
Ansatz c_k = C·q^k/Γ(a k + b), Streuung der q-Schätzer über k; a = 1/κ₁ fest, reelle a, freies komplexes a."""
import numpy as np, mpmath as mp
c = np.load('S1000.npy')
L = mp.mpc('0.31813150520476413', '1.3372357014306895')
mu = mp.sqrt(L); kappa = 2j*mp.pi/(2j*mp.pi - mp.log(mu)); a = 1/kappa
print('a =', a)
K = len(c)
def ratios(b, aa=a):
    out = []
    for k in range(1, K):
        r = (mp.mpc(c[k]) / mp.mpc(c[k-1])) * mp.gamma(aa*k + b) / mp.gamma(aa*(k-1) + b)
        out.append(complex(r))
    return np.array(out)
def spread(params, ks=range(4, 39)):
    b = mp.mpc(params[0], params[1])
    r = ratios(b)
    sel = np.array([r[k-1] for k in ks])
    return np.std(sel) / abs(np.mean(sel)), np.mean(sel)
from scipy.optimize import minimize
best = None
for b0r in [0.5, 1.0, 1.5, 2.0, 3.0]:
    for b0i in [-0.5, 0.0, 0.5]:
        res = minimize(lambda p: spread(p)[0], [b0r, b0i], method='Nelder-Mead', options={'xatol': 1e-8, 'fatol': 1e-12, 'maxiter': 2000})
        if best is None or res.fun < best.fun: best = res
s, m = spread(best.x)
print('Mittag-Leffler-Ansatz c_k = C·q^k/Γ(a k + b), a = 1/κ₁ fest: bestes b = %.4f%+.4fi, rel. Streuung der q-Schätzer (k=4..38) = %.2e, q = %s' % (best.x[0], best.x[1], s, m))
r = ratios(mp.mpc(best.x[0], best.x[1]))
print('q-Schätzer je k:', ' '.join('%d:%.4f%+.4fi' % (k, r[k-1].real, r[k-1].imag) for k in [2, 4, 8, 12, 16, 20, 25, 30, 35, 39]))
# Vergleich: Exponentialtyp-Ansatz a = 1 (c_k = C q^k / Γ(k + b)) und freies reelles a
for aa_test in [1.0, 0.8936, 0.85, 0.95]:
    def sp(p, aa_t=aa_test):
        b = mp.mpc(p[0], p[1]); rr = ratios(b, mp.mpf(aa_t)); sel = rr[3:38]
        return np.std(sel)/abs(np.mean(sel))
    bb = min((minimize(sp, [x, 0.0], method='Nelder-Mead').fun, x) for x in [0.5, 1.0, 2.0])
    print('  reelles a = %.4f: beste rel. Streuung %.2e' % (aa_test, bb[0]))
# freies komplexes a
def spa(p):
    aa = mp.mpc(p[0], p[1]); b = mp.mpc(p[2], p[3]); rr = ratios(b, aa); sel = rr[3:38]
    return np.std(sel)/abs(np.mean(sel))
rb = min((minimize(spa, [0.9, 0.0, x, 0.0], method='Nelder-Mead', options={'maxiter': 4000, 'xatol':1e-9,'fatol':1e-13}) for x in [0.5, 1.0, 2.0]), key=lambda r: r.fun)
print('freies komplexes a: a = %.4f%+.4fi, b = %.4f%+.4fi, Streuung %.2e   (1/κ₁ = %.4f%+.4fi)' % (rb.x[0], rb.x[1], rb.x[2], rb.x[3], rb.fun, float(a.real), float(a.imag)))

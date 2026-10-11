"""D16 — Induktionsschritt N → N+1 für N = 160 … 249 (Fortsetzung von d01, ohne Kreis-Suprema):
Newton-Vorhersage δ = −r w gegen den wahren Schritt c_N, relativer Fehler, |c_N| in ℓ1_{|L|,s}, s_{N+1}, Spitze t_N."""
import sys, time
import numpy as np
from flint import arb, ctx
import d_kern as K
sys.stdout.reconfigure(encoding="utf-8")
lo, hi, st = [int(x) for x in sys.argv[1:4]]
p = K.loesungen(hi + 1, bits=1200)
L = 1.3745570107437075
t0 = time.time()
for N in range(lo, hi + 1, st):
    ctx.prec = 600 + 3 * N
    a = [arb(x.mid()) for x in p[N]]; a1 = [arb(x.mid()) for x in p[N + 1]]
    ap = a + [arb(0)]
    F, J = K.FJ(ap, N + 1)
    r = F[N + 1]
    w = K.loese(J, [arb(0)] * (N + 1) + [arb(1)])
    k = np.arange(N + 2)
    c = np.array([float((a1[i] - ap[i]).mid()) for i in range(N + 2)])
    dl = np.array([float((-r * w[i]).mid()) for i in range(N + 2)])
    out = []
    for s in (0, 4):
        ww = L ** k * (1.0 + k) ** s
        out.append(f"s={s}: |c|={np.abs(c*ww).sum():.3e} rel={np.abs((dl-c)*ww).sum()/np.abs(c*ww).sum():.1e}")
    sig = 1 / float(w[N + 1].mid())
    print(f"N={N} ({time.time()-t0:.0f}s) r={float(r.mid()):+.3e} s={sig:.4e} t={float(a1[N+1].mid()):+.3e} " + " | ".join(out), flush=True)

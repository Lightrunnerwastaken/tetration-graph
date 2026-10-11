"""D15 — Newton von q_N aus: |δ1| = |J(q)^{-1}F(q)|, |δ2| (zweiter Schritt), effektive Konstante |δ2|/|δ1|² in ℓ1_{|L|,s}."""
import json, sys, time
import numpy as np
from flint import arb, ctx
import d_kern as K
sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
kn = json.load(open(K.GEMEINSAM + "/kneser_halbexp_koeff_1000.json"))["a"]
L = 1.3745570107437075
t0 = time.time()
for N in NS:
    ctx.prec = 500 + 3 * N
    q = [arb(str(x)) for x in kn[:N + 1]]
    F, J = K.FJ(q, N); d1 = K.loese(J, F)
    q1 = [arb((q[i] - d1[i]).mid()) for i in range(N + 1)]
    F1, J1 = K.FJ(q1, N); d2 = K.loese(J1, F1)
    k = np.arange(N + 1)
    v1 = np.array([float(x.mid()) for x in d1]); v2 = np.array([float(x.mid()) for x in d2])
    out = []
    for s in (0, 4):
        w = L ** k * (1.0 + k) ** s
        a1, a2 = np.abs(v1 * w).sum(), np.abs(v2 * w).sum()
        out.append(f"s={s}: |δ1|={a1:.3e} |δ2|={a2:.3e} |δ2|/|δ1|²={a2 / a1**2:.3e}")
    M = 8192
    F2 = np.abs(np.fft.ifft(np.concatenate([v2 * L ** k, np.zeros(M - N - 1)])) * M); i = int(np.argmax(F2))
    print(f"N={N} ({time.time()-t0:.0f}s) " + " | ".join(out) + f" | δ2 auf |x|=|L| maximal bei {360*i/M:.2f}°", flush=True)

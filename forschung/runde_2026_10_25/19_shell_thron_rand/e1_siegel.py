"""E1: Siegel-Linearisierung am Shell-Thron-Rand.

Normalform: v = a(z-L), a = ln b, lambda = e^{2 pi i alpha}, b = exp(lambda e^{-lambda}).
g(v) = lambda (e^v - 1) ist konjugiert zu z -> b^z.  Singulärer Wert v = -lambda  (z = 0).
Linearisierer v(w) = w + sum c_n w^n mit v(lambda w) = g(v(w)):  c_n = R_n/(lambda^{n-1}-1).

Aufruf: python e1_siegel.py alpha_name N
"""
import sys
import numpy as np
import mpmath as mp

ALPHAS = {
    "gold": (3 - 5 ** 0.5) / 2,
    "sqrt2": 2 ** 0.5 - 1,
}


def koeff_np(lam, N, s=1.0):
    """skalierte Koeffizienten c_n s^n (gegen Unterlauf)."""
    c = np.zeros(N + 1, complex)
    e = np.zeros(N + 1, complex)
    c[1] = s
    e[0] = 1.0
    e[1] = s
    k = np.arange(N + 1)
    for n in range(2, N + 1):
        # R_n = (1/n) sum_{k=1}^{n-1} k c_k e_{n-k}
        R = np.dot(k[1:n] * c[1:n], e[n - 1:0:-1]) / n
        c[n] = R / (lam ** (n - 1) - 1)
        e[n] = c[n] + R
    return c


def koeff_mp(lam, N):
    c = [mp.mpc(0)] * (N + 1)
    e = [mp.mpc(0)] * (N + 1)
    c[1] = mp.mpc(1)
    e[0] = mp.mpc(1)
    e[1] = mp.mpc(1)
    for n in range(2, N + 1):
        R = mp.fsum(k * c[k] * e[n - k] for k in range(1, n)) / n
        c[n] = R / (lam ** (n - 1) - 1)
        e[n] = c[n] + R
    return c


def radius(c, n0):
    n = np.arange(len(c))
    lg = np.log(np.abs(c[n0:]) + 1e-300)
    # obere Hülle: laufendes Maximum von log|c_n|/n über Fenster
    est = []
    for m in [len(c) // 4, len(c) // 2, len(c) - 1]:
        w = slice(max(n0, m - len(c) // 8), m + 1)
        nn = n[w]
        est.append(float(np.exp(-np.max(np.log(np.abs(c[w])) / nn))))
    # Regression der oberen Hülle (Maxima in Blöcken)
    B = 50
    xs, ys = [], []
    for s in range(n0, len(c) - B, B):
        j = s + int(np.argmax(np.abs(c[s:s + B])))
        xs.append(j)
        ys.append(np.log(abs(c[j])))
    A = np.vstack([np.array(xs), np.log(xs), np.ones(len(xs))]).T
    sol = np.linalg.lstsq(A, np.array(ys), rcond=None)[0]
    return est, float(np.exp(-sol[0])), float(sol[1])


def birkhoff(lam, p, N, dps):
    """gewichtetes Birkhoff-Mittel von lambda^{-k} g^k(p)  -> phi(p)."""
    with mp.workdps(dps):
        lam = mp.mpc(lam)
        v = mp.mpc(p)
        num = mp.mpc(0)
        den = mp.mpf(0)
        li = 1 / lam
        f = mp.mpc(1)
        maxabs = 0
        for k in range(N):
            t = mp.mpf(k + 0.5) / N
            w = mp.exp(-1 / (t * (1 - t)))
            num += w * f * v
            den += w
            maxabs = max(maxabs, abs(v))
            v = lam * (mp.exp(v) - 1)
            f *= li
        return num / den, maxabs


if __name__ == "__main__":
    name = sys.argv[1] if len(sys.argv) > 1 else "gold"
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
    al = ALPHAS[name]
    with mp.workdps(40):
        alm = (3 - mp.sqrt(5)) / 2 if name == "gold" else mp.sqrt(2) - 1
        lamm = mp.exp(2j * mp.pi * alm)
        bm = mp.exp(lamm * mp.exp(-lamm))
        print("alpha", mp.nstr(alm, 20), " lambda", mp.nstr(lamm, 20), " b", mp.nstr(bm, 20))
    lam = complex(lamm)
    S = 1.7
    c = koeff_np(lam, N, S)
    cm = koeff_mp(lamm, 300)
    rel = max(abs(complex(cm[n]) * S**n - c[n]) / abs(complex(cm[n]) * S**n) for n in range(2, 301))
    print("numpy vs mpmath (n<=300) max rel. Abweichung", rel)
    est, rfit, pw = radius(c, 20)
    est = [S * e for e in est]; rfit = S * rfit
    print("Wurzeltest-Fenster r ~", est, " Hüllen-Fit r =", rfit, " Potenz n^", pw)
    np.save(f"cache/c_{name}_{N}.npy", c)
    for NB in ([16000] if len(sys.argv) < 4 else []):
        for dps in ([30, 60] if NB <= 2000 else [25]):
            w0, mx = birkhoff(lamm, -lamm, NB, dps)
            print(f"Birkhoff N={NB} dps={dps}: phi(-lam) = {mp.nstr(w0, 15)}  |.|={mp.nstr(abs(w0), 15)}  max|orbit|={mp.nstr(mx, 6)}")

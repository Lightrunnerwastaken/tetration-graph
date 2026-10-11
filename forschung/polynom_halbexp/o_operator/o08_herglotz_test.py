"""O08 — Prüfung des Herglotz-Kriteriums (Satz O2) im symmetrischen elliptischen Modell:
K akkretiv auf Polynomen in H²(D) ⟺ Re F ≥ 0 auf |Y|<1, F(Y) = (1−2tY)²/(1−tY), t = sin(c/2) ⟺ c ≤ c* = 2 arcsin(√((4√2−5)/4)).
Hochgenau: kleinster Eigenwert des hermiteschen Teils H_N = (K_N + K_N*)/2 (in mpmath), dazu min Re F auf dem Rand."""
import sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 120
cstar = 2 * mp.asin(mp.sqrt((4 * mp.sqrt(2) - 5) / 4))
print("c* =", mp.nstr(cstar, 15))

def H_matrix(c, N):
    tau = mp.tan(mp.mpf(c) / 2)
    # f = (x+τ)/(1−τx); a = 1/m = (1−τ²−2τx)²/((1+τ²)(1−τx)²)
    n = N + 1
    geo = [tau ** k for k in range(n)]                      # 1/(1−τx) = Σ τ^k x^k
    fk = [tau] + [geo[k] + tau * geo[k - 1] * tau for k in range(1, n)]   # (x+τ)Σ τ^k x^k: [x^k] = τ^{k−1} + τ·τ^k
    fk = [tau * geo[0]] + [geo[k - 1] + tau * geo[k] for k in range(1, n)]
    def mul(p, q):
        r = [mp.mpf(0)] * n
        for i in range(n):
            if p[i] == 0: continue
            for j in range(n - i):
                r[i + j] += p[i] * q[j]
        return r
    num = [1 - tau ** 2, -2 * tau] + [mp.mpf(0)] * (n - 2)
    num2 = mul(num, num)
    g2 = mul(geo, geo)
    a = [x / (1 + tau ** 2) for x in mul(num2, g2)]
    pw = [[mp.mpf(1)] + [mp.mpf(0)] * N]
    for _ in range(N):
        pw.append(mul(pw[-1], fk))
    Km = mp.matrix(n, n)
    for j in range(n):
        col = mul(pw[j], a)
        for k in range(n):
            Km[k, j] = col[k]
    H = (Km + Km.T) / 2
    return H

for c in (0.80, 0.83, 0.84, 0.86):
    t = mp.sin(mp.mpf(c) / 2)
    minRe = min(mp.re((1 - 2 * t * mp.expj(p)) ** 2 / (1 - t * mp.expj(p))) for p in mp.linspace(0, mp.pi, 2001))
    res = []
    for N in (10, 20, 30, 40):
        ev = mp.eigsy(H_matrix(c, N), eigvals_only=True)
        res.append(mp.nstr(min(ev), 4))
    print(f"c={c}: min Re F auf |Y|=1: {mp.nstr(minRe, 5)};  min eig H_N (N=10,20,30,40): {res}", flush=True)

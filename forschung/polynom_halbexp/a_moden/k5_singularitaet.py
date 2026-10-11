"""Strang A, Schwerpunkt (Leitung 20:00): Knesers Singularitätsentwicklung der halben e-Funktion an L.

Teil A (exakt/formal): Schröder-Funktion an L. φ = Ψ^{-1}: φ(Ls) = exp(φ(s)), φ(0) = L, φ'(0) = 1 (ganze Funktion).
    Reguläre halbe Iterierte f_reg(z) = φ(√L·Ψ(z)), √L = e^{L/2}; ihre Taylor-Koeffizienten an L: c_m.
Teil C (numerisch, Motor): Knesers θ-Abbildung an L. Mit u-Koordinate s = e^{Lu}, z = φ(s):
    A(z) = u + ϑ(u),  ϑ(u) = Σ_{m≥0} κ_m e^{2πimu}  (1-periodisch, beschränkt für Im u → +∞).
    κ_m aus Fourier-Mitteln über u ∈ [u0, u0+1) + iv (zwei Höhen v, negative m müssen verschwinden).
Folgerung (Rechnung im Bericht): f_K(z) = φ(√L·s·e^{L ε(q)}), q = s^γ = e^{2πiu}, γ = 2πi/L, mit
    ε(q) = ϑ(u) − ϑ(U), U = u + ½ + ε  ⇒  ε = 2κ_1 q + 4πiκ_1² q² + (2κ_3 − 8πiκ_1κ_2 − 16π²κ_1³) q³ + …
    also f_K − f_reg = 2 L^{3/2} κ_1 s^{1+γ}(1 + O(s)) + O(s^{1+2γ}).
Prüfung: f_K(z) (Motor) gegen diese Formel an Punkten nahe L.
"""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402
from fatou_backend import FatouGP  # noqa: E402

HIER = Path(__file__).resolve().parent
GEM = HIER.parent / "gemeinsam"
DATEN = HIER / "daten"
mp.mp.dps = 70
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
SL = mp.exp(L / 2)                      # √L (Hauptzweig, = f_K'(L))
G = 2j * mp.pi / L                      # γ
K = 90                                  # Reihenlänge


# ---------- Teil A: φ = Ψ^{-1}
def exp_reihe(g, n):
    """exp(g(s)) für g mit g_0 = 0, als Reihe bis s^{n-1}."""
    E = [mp.mpc(1)] + [mp.mpc(0)] * (n - 1)
    for k in range(1, n):
        E[k] = mp.fsum(j * g[j] * E[k - j] for j in range(1, k + 1)) / k
    return E


phi = [L, mp.mpc(1)] + [mp.mpc(0)] * (K - 2)
E = [mp.mpc(1), mp.mpc(1)] + [mp.mpc(0)] * (K - 2)
for k in range(2, K):
    R = mp.fsum(j * phi[j] * E[k - j] for j in range(1, k)) / k     # ohne den j = k-Term (phi_k unbekannt)
    phi[k] = L * R / (L ** k - L)
    E[k] = phi[k] + R
# Prüfung φ(Ls) = exp(φ(s)) bei s = 0.3 e^{i}
s0 = mp.mpf("0.3") * mp.expj(1)
pv = lambda c, x: mp.fsum(c[k] * x ** k for k in range(len(c)))  # noqa: E731
print(f"Teil A: |φ(Ls) − exp(φ(s))| bei |s| = 0.3: {mp.nstr(abs(pv(phi, L * s0) - mp.exp(pv(phi, s0))), 3)}; "
      f"|φ_k| bei k = 10, 30, 60, 89: {[mp.nstr(abs(phi[k]), 3) for k in (10, 30, 60, 89)]}")


def mul(a, b, n):
    return [mp.fsum(a[i] * b[k - i] for i in range(k + 1)) for k in range(n)]


# Ψ = φ^{-1} als Reihe in t = z − L: Ψ(t) = t + …, über Newton-Iteration  Ψ ← Ψ − (φ(Ψ) − L − t)/φ'(Ψ)
NA = 50
psi = [mp.mpc(0), mp.mpc(1)] + [mp.mpc(0)] * (NA - 2)


def komp(c, inner, n):
    """c(inner(t)) als Reihe bis t^{n-1}, inner_0 = 0."""
    out = [mp.mpc(0)] * n
    pot = [mp.mpc(1)] + [mp.mpc(0)] * (n - 1)
    for k in range(min(len(c), n)):
        for i in range(n):
            out[i] += c[k] * pot[i]
        pot = mul(pot, inner, n)
    return out


for it in range(7):
    ph = komp(phi, psi, NA)
    dphi = [k * phi[k] for k in range(1, len(phi))]
    dph = komp(dphi, psi, NA)
    res = list(ph)
    res[0] -= L
    res[1] -= 1
    # Division res/dph als Reihe
    q = [mp.mpc(0)] * NA
    for k in range(NA):
        q[k] = (res[k] - mp.fsum(q[i] * dph[k - i] for i in range(k))) / dph[0]
    psi = [a - b for a, b in zip(psi, q)]
    print(f"   Ψ-Newton {it}: max|Rest| = {mp.nstr(max(abs(x) for x in res), 3)}")
    if max(abs(x) for x in res) < mp.mpf(10) ** -60:
        break
freg = komp([phi[k] * SL ** k for k in range(len(phi))], psi, NA)   # f_reg(L + t) = Σ freg_m t^m
print(f"f_reg an L: f(L) = {mp.nstr(freg[0], 15)}, f'(L) = {mp.nstr(freg[1], 15)} (√L = {mp.nstr(SL, 15)})")
cm = [freg[m + 1] / SL for m in range(NA - 1)]    # f_reg − L = √L t (1 + Σ_{m≥1} c_m t^m)
print("c_m (f_reg − L = √L t(1 + Σ c_m t^m)): " + ", ".join(f"c_{m} = {mp.nstr(cm[m], 12)}" for m in range(1, 6)))
print(f"   |c_m|^(1/m) bei m = 10, 20, 40: {[mp.nstr(abs(cm[m]) ** (mp.mpf(1) / m), 5) for m in (10, 20, 40)]}")

# ---------- Teil C: κ_m aus dem Motor
gp = FatouGP(dps=57, fatou_gp="fork")
s_ = gp.session("exp(1)")
m = 64
kap = {}
for v in (mp.mpf(1), mp.mpf("1.5")):
    # Bogen auf der Seite von 0: arg s = Im(L u) ≈ −1.8 in der Mitte
    u0 = (mp.mpf("-1.8") - mp.im(L) * 0 - mp.re(L) * v) / mp.im(L) - mp.mpf("0.5")
    us = [mp.mpc(u0 + mp.mpf(j) / m, v) for j in range(m)]
    ss = [mp.exp(L * u) for u in us]
    zs = [pv(phi, s) for s in ss]
    A = s_.slog_batch(zs)
    th = [a - u for a, u in zip(A, us)]
    F = {k: mp.fsum(t * mp.exp(-2j * mp.pi * k * u) for t, u in zip(th, us)) / m for k in range(-3, 7)}
    kap[v] = F
    print(f"\nTeil C, Höhe v = {v}: |s| ∈ [{mp.nstr(min(abs(s) for s in ss), 3)}, {mp.nstr(max(abs(s) for s in ss), 3)}], "
          f"max|z − L| = {mp.nstr(max(abs(z - L) for z in zs), 3)}")
    for k in range(-3, 7):
        print(f"   κ_{k:+d} = {mp.nstr(F[k], 25)}")
v1, v2 = list(kap)
print("\nHöhenvergleich |κ_m(v=1) − κ_m(v=1.5)|/|κ_m|: " + ", ".join(
    f"m={k}: {mp.nstr(abs(kap[v1][k] - kap[v2][k]) / abs(kap[v2][k]), 2)}" for k in range(0, 5)))
kp = kap[v2]

# ε(q) als Reihe bis q^4 (Fixpunkt), dann Prüfung f_K = φ(√L s e^{Lε})
def eps_reihe(kp, n=6):
    e = [mp.mpc(0)] * n
    for _ in range(n + 2):
        # ε = Σ_m κ_m q^m (1 − (−1)^m e^{2πimε})
        neu = [mp.mpc(0)] * n
        for mm in range(1, n):
            ex = exp_reihe([2j * mp.pi * mm * x for x in e], n)       # e^{2πimε(q)}
            fak = [(1 if k == 0 else 0) - (-1) ** mm * ex[k] for k in range(n)]
            for k in range(n - mm):
                neu[k + mm] += kp[mm] * fak[k] if mm in kp else 0
        e = neu
    return e


eps = eps_reihe(kp)
print("\nε(q) = " + " + ".join(f"({mp.nstr(eps[k], 12)}) q^{k}" for k in range(1, 5)))
print(f"   Soll 2κ_1 = {mp.nstr(2 * kp[1], 12)}, 4πiκ_1² = {mp.nstr(4j * mp.pi * kp[1] ** 2, 12)}")
# Prüfpunkte nahe L (eigene u, Höhe 1.2)
us = [mp.mpc(mp.mpf(x), mp.mpf("1.2")) for x in ("-1.9", "-1.6", "-1.3", "-1.0")]
ss = [mp.exp(L * u) for u in us]
zs = [pv(phi, s) for s in ss]
fK = s_.sexp_batch([a + mp.mpf("0.5") for a in s_.slog_batch(zs)])
gp.close()
print("\nPrüfung f_K(z) (Motor) gegen φ(√L s e^{Lε(q)}) und gegen f_reg + 2L^{3/2}κ_1 s^{1+γ}:")
for u, s, z, fk in zip(us, ss, zs, fK):
    q = mp.exp(2j * mp.pi * u)
    e = mp.fsum(eps[k] * q ** k for k in range(1, len(eps)))
    voll = pv(phi, SL * s * mp.exp(L * e))
    fr = pv(phi, SL * s)
    lead = fr + 2 * L * SL * kp[1] * s * q
    print(f"   |z−L| = {mp.nstr(abs(z - L), 3)}: |f_K − f_reg| = {mp.nstr(abs(fk - fr), 3)}, |f_K − Formel(ε bis q^5)| = "
          f"{mp.nstr(abs(fk - voll), 3)}, |f_K − (f_reg + Leitterm)| = {mp.nstr(abs(fk - lead), 3)}")
json.dump({"kappa": {str(k): [mp.nstr(mp.re(x), 40), mp.nstr(mp.im(x), 40)] for k, x in kp.items()},
           "c_m": [[mp.nstr(mp.re(x), 40), mp.nstr(mp.im(x), 40)] for x in cm[:40]],
           "phi": [[mp.nstr(mp.re(x), 40), mp.nstr(mp.im(x), 40)] for x in phi[:60]],
           "psi": [[mp.nstr(mp.re(x), 40), mp.nstr(mp.im(x), 40)] for x in psi[:40]],
           "eps_q": [[mp.nstr(mp.re(x), 30), mp.nstr(mp.im(x), 30)] for x in eps]},
          open(DATEN / "singularitaet_theta.json", "w", encoding="utf-8"), indent=0)
print("geschrieben", time.strftime("%H:%M"))

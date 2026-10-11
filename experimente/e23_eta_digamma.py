"""exp.eta_digamma — Basis η: Kehrwert-Koordinate, Konstante C_T, Bruchsumme mit Digamma; Kehrwert-Turmsummen über η.

 q = (e − T)/(2e):  q(x+1) = (1 − e^{−2q(x)})/2;   Z = 1/q:  Z(x+1) = 1 + coth(1/Z(x)).
 Fatou-Koordinate Φ(Z) = Z − (1/3) log Z + Σ c_k Z^{−k} (c_k rational), Φ(Z(x+1)) = Φ(Z(x)) + 1.
 C_T := lim Φ(Z(n)) − n  (nur ganze Höhen, T(0) = 1)  ⇒  Z(x) = x + (1/3) log x + C_T + o(1).
 Bruchsumme S(x) = e x + Σ_{k≥1} [T(k) − T(k+x)] = e x − 2e[ψ(x+1) + γ] − 2e R_T(x),
 R_T(∞) = Σ_k [q(k) − 1/k]  ⇒  S(x) = e x − 2e log x − 2eγ + Σ_k [T(k) − e + 2e/k] + o(1).
 Über η: C_b = Σ 1/T_b(k), R_b(x) = Σ [1/T_b(k) − 1/T_b(k+x)], C_b − R_b(x) = Σ_{m≥0} b^{−T_b(x+m)}.
"""
import mpmath as mp
import sympy as sp

from gemeinsam import Lauf, kneser, plt

mp.mp.dps = 60
lauf = Lauf("exp.eta_digamma", __file__)
E = mp.e
ETA = mp.e ** (1 / mp.e)

# ------------------------------------------------------------------ Fatou-Reihe (exakt, rational)
K = 24
t = sp.symbols("t")
eps = sp.series(sp.coth(t) - 1 / t, t, 0, K + 3).removeO()
cs = sp.symbols(f"c1:{K + 1}")
G = eps - sp.Rational(1, 3) * sp.log(1 + t + t * eps)
for k in range(1, K + 1):
    G += cs[k - 1] * t ** k * ((1 + t + t * eps) ** (-k) - 1)
G = sp.expand(sp.series(G, t, 0, K + 2).removeO())
loes = {}
for k in range(1, K + 1):
    gl = sp.expand(G.coeff(t, k + 1).subs(loes))
    loes[cs[k - 1]] = sp.solve(gl, cs[k - 1])[0]
c_fatou = [mp.mpf(sp.Rational(loes[cs[k - 1]]).p) / mp.mpf(sp.Rational(loes[cs[k - 1]]).q) for k in range(1, K + 1)]
print("  c_k:", [str(loes[cs[k]]) for k in range(4)])


def Phi(Z):
    return Z - mp.log(Z) / 3 + mp.fsum(c_fatou[k - 1] * Z ** (-k) for k in range(1, K + 1))


def Phi_strich(Z):
    return 1 - 1 / (3 * Z) - mp.fsum(k * c_fatou[k - 1] * Z ** (-k - 1) for k in range(1, K + 1))


def Z_von(y):
    return 2 * E / (E - y)


def f(y):
    return mp.exp(y / E)


# ------------------------------------------------------------------ ganze Höhen, C_T
NMAX = 20000
bahn = [mp.mpf(1)]
for _ in range(NMAX):
    bahn.append(f(bahn[-1]))
q = [(E - y) / (2 * E) for y in bahn]
res = max(max(abs(q[n + 1] - (1 - mp.exp(-2 * q[n])) / 2), abs(1 / q[n + 1] - 1 - mp.coth(q[n]))) for n in (0, 5, 100, 5000))
lauf.pruefe("exakte Umformungen: q(x+1) = (1 − e^{−2q})/2 und Z(x+1) = 1 + coth(1/Z)", res < mp.mpf(10) ** -50, mp.nstr(res, 3))
Ct = {n: Phi(Z_von(bahn[n])) - n for n in (2500, 5000, 10000, 20000)}
lauf.pruefe("C_T = lim Φ(Z(n)) − n ist stabil (n = 2500 … 20000)", abs(Ct[2500] - Ct[20000]) < mp.mpf(10) ** -35,
            f"C_T = {mp.nstr(Ct[20000], 30)}")
C_T = Ct[20000]
abw = {n: Z_von(bahn[n]) - n - mp.log(n) / 3 - C_T for n in (100, 1000, 10000)}
lauf.pruefe("Asymptotik Z(n) = n + (1/3) log n + C_T + O(log n / n)", all(abs(abw[n]) < mp.log(n) / n for n in abw) and abs(abw[10000]) < abs(abw[100]),
            {n: mp.nstr(v, 4) for n, v in abw.items()})


def Z_aus_hoehe(x):
    """Z(x) aus Φ(Z) = x + C_T (Newton), gültig für grosse x."""
    Z = x + mp.log(x) / 3 + C_T
    for _ in range(50):
        dZ = (Phi(Z) - x - C_T) / Phi_strich(Z)
        Z -= dZ
        if abs(dZ) < mp.mpf(10) ** (-mp.mp.dps + 5):
            break
    return Z


# ------------------------------------------------------------------ gebrochene Höhen über die Fatou-Koordinate
def alpha(y, N=3000):
    for _ in range(N):
        y = f(y)
    return Phi(Z_von(y)) - N


def T_bruch(x0):
    """T(x0) für x0 in [0, 1): α(y) = x0 + C_T."""
    return mp.findroot(lambda y: alpha(y) - x0 - C_T, (mp.mpf(1), ETA), solver="anderson")


mp.mp.dps = 50
xh = mp.mpf("0.5")
Th = T_bruch(xh)
kontrolle = abs(alpha(Th, N=1000) - alpha(Th, N=3000))
lauf.pruefe("Fatou-Koordinate unabhängig von der Iterationstiefe (N = 1000 vs 3000)", kontrolle < mp.mpf(10) ** -30, mp.nstr(kontrolle, 3))


def bruchsumme(x_ganz, x_frac, Ty, KS=20000):
    """S(x) = e x + Σ_{k≥1}[T(k) − T(k+x)], x = x_ganz + x_frac, Schwanz nach Euler–Maclaurin."""
    x = x_ganz + x_frac
    y = Ty
    for _ in range(x_ganz):
        y = f(y)
    s = mp.mpf(0)
    for k in range(1, KS + 1):
        y = f(y)
        s += bahn[k] - y
    qa = lambda m: 1 / Z_aus_hoehe(m)
    g = lambda m: 2 * E * (qa(m + x) - qa(m))                         # T(m) − T(m+x)
    integral = 2 * E * mp.quad(qa, [KS, KS + x])                     # ∫_K^∞ g = 2e ∫_K^{K+x} q
    gstrich = (g(KS + mp.mpf("0.5")) - g(KS - mp.mpf("0.5")))
    schwanz = -integral - g(KS) / 2 - gstrich / 12
    return E * x + s + schwanz


S_h = bruchsumme(0, xh, Th)
S_h1 = bruchsumme(1, xh, Th)
lauf.pruefe("Bruchsumme: S(3/2) − S(1/2) = T(3/2)", abs(S_h1 - S_h - f(Th)) < mp.mpf(10) ** -14, mp.nstr(abs(S_h1 - S_h - f(Th)), 3))
S3 = bruchsumme(3, mp.mpf(0), mp.mpf(1))
lauf.pruefe("Bruchsumme an ganzer Stelle: S(3) = T(1) + T(2) + T(3)", abs(S3 - bahn[1] - bahn[2] - bahn[3]) < mp.mpf(10) ** -14,
            mp.nstr(abs(S3 - bahn[1] - bahn[2] - bahn[3]), 3))

# Konstante der Asymptotik aus ganzen Höhen
KN = NMAX
reg = mp.fsum(bahn[k] - E + 2 * E / k for k in range(1, KN + 1))
qa = lambda m: 1 / Z_aus_hoehe(m)
h_reg = lambda m: -2 * E * qa(m) + 2 * E / m                          # T(m) − e + 2e/m
schw = mp.quad(h_reg, [KN, mp.inf]) - h_reg(KN) / 2 - (h_reg(KN + mp.mpf("0.5")) - h_reg(KN - mp.mpf("0.5"))) / 12
reg_wert = reg + schw
reg_10k = mp.fsum(bahn[k] - E + 2 * E / k for k in range(1, 10001))
schw_10k = mp.quad(h_reg, [10000, mp.inf]) - h_reg(10000) / 2 - (h_reg(mp.mpf("10000.5")) - h_reg(mp.mpf("9999.5"))) / 12
lauf.pruefe("regularisierte Reihe Σ_k [T(k) − e + 2e/k] mit Schwanz stabil (10⁴ vs 2·10⁴ Glieder)", abs(reg_10k + schw_10k - reg_wert) < mp.mpf(10) ** -14,
            f"Wert {mp.nstr(reg_wert, 20)}")
K_eta = -2 * E * mp.euler + reg_wert
dig = {}
for xg in (20, 200):
    Sx = bruchsumme(xg, xh, Th)
    xx = xg + xh
    D = Sx - (E * xx - 2 * E * (mp.digamma(xx + 1) + mp.euler)) - reg_wert          # = −2e (R_T(x) − R_∞)
    vorhersage = -2 * E * ((mp.log(xx) + 1) / 3 + C_T) / xx                         # aus q(m) ≈ 1/m − ((1/3) log m + C_T)/m²
    dig[xg] = (D, D / vorhersage - 1)
lauf.pruefe("Digamma-Zerlegung: S(x) = e x − 2e(ψ(x+1)+γ) + Σ_k[T(k) − e + 2e/k] − 2e((log x + 1)/3 + C_T)/x + …",
            abs(dig[200][1]) < 0.05 and abs(dig[200][1]) < abs(dig[20][1]), {k: f"Rest {mp.nstr(v[0], 4)}, rel. zur Vorhersage {mp.nstr(v[1], 3)}" for k, v in dig.items()})

# Ableitungsreihe Σ T'(k): T'(0) = 1/α'(1), T'(k) = T'(0) Π_{j≤k} T(j)/e
mp.mp.dps = 60
y, prod = mp.mpf(1), mp.mpf(1)
for _ in range(3000):
    y_neu = f(y)
    prod *= y_neu / E
    y = y_neu
alpha_strich_1 = Phi_strich(Z_von(y)) * 2 * E / (E - y) ** 2 * prod
T0s = 1 / alpha_strich_1
Ts = [T0s]
for k in range(1, NMAX + 1):
    Ts.append(Ts[-1] * bahn[k] / E)
asym = lambda m: 2 * E / (Z_aus_hoehe(m) ** 2 * Phi_strich(Z_aus_hoehe(m)))
lauf.pruefe("T′(n) über die Transferformel = T′ aus der Fatou-Asymptotik (n = 20000)", abs(Ts[NMAX] / asym(mp.mpf(NMAX)) - 1) < mp.mpf(10) ** -25,
            mp.nstr(abs(Ts[NMAX] / asym(mp.mpf(NMAX)) - 1), 3))
summe_abl = mp.fsum(Ts[1:]) + mp.quad(asym, [NMAX, mp.inf]) - asym(mp.mpf(NMAX)) / 2
lauf.pruefe("Σ_{k≥1} T′(k) konvergiert (T′(k) ≈ 2e/k²)", abs(Ts[1000] * 1000 ** 2 / (2 * E) - 1) < 0.02, f"Σ T′(k) = {mp.nstr(summe_abl, 18)}")

# ------------------------------------------------------------------ Kehrwert-Turmsummen über η
mp.mp.dps = 40
Cb = {}
for name, bb in (("2", mp.mpf(2)), ("e", mp.e), ("10", mp.mpf(10))):
    y, s = mp.mpf(1), mp.mpf(0)
    for _ in range(6):
        y = mp.power(bb, y) if y < 10 ** 6 else mp.inf
        s += 1 / y
    Cb[name] = s
lauf.pruefe("C_2 = 1/2 + 1/4 + 1/16 + 1/65536 + … = 0.8125152587890625 + 2^{−65536} + …", abs(Cb["2"] - mp.mpf("0.8125152587890625")) < mp.mpf(10) ** -38,
            mp.nstr(Cb["2"], 20))
T2h = mp.re(kneser("sexp", "2", [mp.mpf("0.5")])[0])
werte = {}
for xg in (0, 1, 2):
    y = T2h
    for _ in range(xg):
        y = mp.power(2, y)
    rest, z = mp.mpf(0), y
    for _ in range(5):
        z = mp.power(2, z) if z < 10 ** 6 else mp.inf
        rest += 1 / z
    werte[xg] = rest / mp.power(2, -y)                               # (C − R(x)) / b^{−T(x)}
lauf.pruefe("über η: C_b − R_b(x) = Σ_m b^{−T(x+m)}, erster Term dominiert (Verhältnis → 1, Basis 2, x = 1/2, 3/2, 5/2)",
            abs(werte[2] - 1) < abs(werte[0] - 1) and abs(werte[2] - 1) < 1e-20, {k: mp.nstr(v, 6) for k, v in werte.items()})

lauf.daten = {"c_fatou": [str(loes[cs[k]]) for k in range(6)], "C_T": mp.nstr(C_T, 35), "T_eta_halb": mp.nstr(Th, 30),
              "S_eta_halb": mp.nstr(S_h, 18), "regularisiert": mp.nstr(reg_wert, 20), "K_eta": mp.nstr(K_eta, 20),
              "T_strich_0": mp.nstr(T0s, 25), "summe_T_strich": mp.nstr(summe_abl, 18),
              "C_b": {k: mp.nstr(v, 25) for k, v in Cb.items()}, "T_2_halb": mp.nstr(T2h, 25)}

fig, ax = plt.subplots(1, 2, figsize=(10.5, 4))
ns = [10 * 2 ** k for k in range(0, 11)]
ax[0].semilogx(ns, [float(Z_von(bahn[n]) - n - mp.log(n) / 3) for n in ns], "o-", ms=3, label="Z(n) − n − (1/3) log n")
ax[0].axhline(float(C_T), ls="--", color="#b5523b", label=f"C_T = {mp.nstr(C_T, 10)}")
ax[0].set_xlabel("n"); ax[0].set_title("Basis η: Kehrwert-Koordinate"); ax[0].legend(fontsize=8)
xs = [mp.mpf(i) / 4 for i in range(2, 41)]
Tx = {}
for xv in xs:
    ganz, frac = int(mp.floor(xv)), xv - int(mp.floor(xv))
    y = Th if frac == mp.mpf("0.5") else (mp.mpf(1) if frac == 0 else None)
    if y is None:
        continue
    for _ in range(ganz):
        y = f(y)
    Tx[xv] = y
ax[1].plot([float(x) for x in Tx], [float(E - v) for v in Tx.values()], "o", ms=3, label="e − T_η(x) (ganz und halbzahlig)")
ax[1].plot([float(x) for x in xs], [float(2 * E / (x + mp.log(x) / 3 + C_T)) for x in xs], "-", lw=1, label="2e/(x + (1/3)log x + C_T)")
ax[1].set_xlabel("x"); ax[1].set_title("Annäherung an e"); ax[1].legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = (f"Basis η: C_T = {mp.nstr(C_T, 20)} (aus ganzen Höhen), T_η(1/2) = {mp.nstr(Th, 16)}, "
                        f"Σ_k[T(k) − e + 2e/k] = {mp.nstr(reg_wert, 16)}, Σ T′(k) = {mp.nstr(summe_abl, 14)}; Digamma-Zerlegung bestätigt. "
                        f"C_2 = {mp.nstr(Cb['2'], 17)}, C_e = {mp.nstr(Cb['e'], 17)}.")
lauf.schreibe()

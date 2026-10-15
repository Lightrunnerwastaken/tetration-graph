"""Lücken S3/S4/Z0 (Basis e): Zahlen zu BERICHT_LUECKEN_ANALYTISCH.md, Abschnitte 2–4.

Teil A  Taylor von sexp bei 0 (etappe3, 130 Stellen): Univalenz auf |z| < r (Noshiro–Warschawski),
        Schranken e0, e1 für sexp z = 1 + A z(1+ε0), log sexp z = A z(1+ε1) auf |z| ≤ r.
Teil B  Ein-Punkt-Zertifikat für die Zweigwahl (S4): P = e^{L(z0+θ(z0))} liegt in σ(ℍ), d. h. σ(Ψ(P)) = P,
        mit Trennradius ρ = |σ'(w0)| Im w0 / 4 (Koebe-1/4 auf B(w0, Im w0) ⊂ ℍ).
Teil C  Probe der Grunsky-Ungleichungen |log Q(w)| ≤ d(w), |Λ(w) + log(1−|ζ|²)| ≤ d(w).
Teil D  Konstanten der Majoranten auf D_r und Vergleich mit θ selbst (θ = 3 − z + h(u3) + 2πik/L).
Kein Motor. Aufruf: python l1_schranken.py   (< 1 min)
"""
import json
import sys
from pathlib import Path

import mpmath as mp

HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER.parent / "etappe3"))
from lib3 import Basis, ThetaReihe, lade_theta  # noqa: E402

mp.mp.dps = 50
B = Basis(K=60)
L = B.L
Lb = mp.conj(L)
ImL = mp.im(L)
out = {}

# ---------------- Teil A
tay = json.loads((HIER.parent / "etappe3" / "daten" / "taylor_130_256.json").read_text())
c = [mp.mpf(x) if "j" not in x else mp.mpc(complex(x)) for x in tay["zentren"]["0.0"]]
c = [mp.mpc(v) for v in c]
A = mp.re(c[1])
N = len(c) - 1


def sexp0(z, abl=False):
    s = mp.fsum(c[k] * z ** k for k in range(N + 1))
    if not abl:
        return s
    return s, mp.fsum(k * c[k] * z ** (k - 1) for k in range(1, N + 1))


r = mp.mpf("0.1")
NW = mp.fsum(k * abs(c[k]) * r ** (k - 1) for k in range(2, N + 1))
# grobe Restschranke: |c_k| ≤ M1 (Kreis |z| = 1, Taylor-Polynom), Rest Σ_{k>N} k M1 r^{k-1}
M1 = max(abs(sexp0(mp.expjpi(2 * mp.mpf(j) / 256))) for j in range(256))
rest = mp.fsum(k * M1 * r ** (k - 1) for k in range(N + 1, N + 200))
print(f"A = {mp.nstr(A, 20)};  Σ_(k≥2) k|c_k| r^(k-1) = {mp.nstr(NW, 8)} (+Rest ≤ {mp.nstr(rest, 3)}) < A  ⇒ sexp univalent auf |z| < {r}")
e0 = e1 = es = mp.mpf(0)
for j in range(720):
    z = r * mp.expjpi(2 * mp.mpf(j) / 720)
    s, sd = sexp0(z, True)
    e0 = max(e0, abs((s - 1) / (A * z) - 1))
    e1 = max(e1, abs(mp.log(s) / (A * z) - 1))
    es = max(es, abs(sd / A - 1))
print(f"auf |z| ≤ r: |ε0| ≤ {mp.nstr(e0, 5)}, |ε1| ≤ {mp.nstr(e1, 5)}, |sexp'/A − 1| ≤ {mp.nstr(es, 5)}")
out["A"] = {"r": str(r), "NW": mp.nstr(NW, 10), "e0": mp.nstr(e0, 6), "e1": mp.nstr(e1, 6), "es": mp.nstr(es, 6)}

# ---------------- Teil B: Zweig-Zertifikat
a = lade_theta(nmax=4096)
T = ThetaReihe(a, B)
zert = []
for z0 in [mp.mpc(0, "0.05"), mp.mpf("0.07") * mp.expjpi(mp.mpf(1) / 3), mp.mpf("0.08") * mp.expjpi(mp.mpf(3) / 4),
           mp.mpc("0.5", "0.05")]:
    th, th1 = T.theta(z0, 1)
    P = mp.exp(L * (z0 + th))
    w0, dPsi = B.Psi(P, True)
    sw, n, su, dlog = B.sigma(w0, True)
    rho = abs(dlog * sw) * mp.im(w0) / 4
    diff = abs(sw - P)
    st = abs(w0 - sexp0(z0)) if abs(z0) < r else mp.nan
    # Zweigzahl k: θ(z0) − (3 − z0 + h(u3)) = 2πik/L
    u = w0
    for _ in range(3):
        u = mp.log(u)
    k = (th - (3 - z0 + B.h(u))) * L / (2j * mp.pi)
    zert.append({"z0": mp.nstr(z0, 6), "|sigma(Psi(P))-P|": mp.nstr(diff, 3), "rho": mp.nstr(rho, 4),
                 "|Psi(P)-sexp_Taylor|": mp.nstr(st, 3), "k": mp.nstr(k, 8), "Im w0": mp.nstr(mp.im(w0), 4)})
    print(f"z0 = {mp.nstr(z0, 6)}:  |σ(Ψ(P)) − P| = {mp.nstr(diff, 3)}  <  ρ = {mp.nstr(rho, 4)};"
          f"  |Ψ(P) − sexp_Taylor(z0)| = {mp.nstr(st, 3)};  Zweigzahl k = {mp.nstr(k, 10)}")
out["B"] = zert
z05 = mp.mpc(0, "0.05")
k1 = mp.nint(mp.re((T.theta(z05) - (3 - z05 + B.h(mp.log(mp.log(mp.log(sexp0(z05))))))) * L / (2j * mp.pi)))
print("Zweigzahl k1 =", k1)


# ---------------- Teil C: Grunsky-Probe
def zeta(w):
    return (w - L) / (w - Lb)


def dhyp(w):
    return mp.acosh(1 + abs(w - L) ** 2 / (2 * mp.im(w) * ImL))


def stetig_log(f, w, schritte=400):
    """stetiger Zweig von log f längs [L+ε(w−L), w], Start beim Hauptzweig."""
    prev = None
    for j in range(1, schritte + 1):
        t = mp.mpf(j) / schritte
        v = mp.log(f(L + t * (w - L)))
        if prev is not None:
            v += 2j * mp.pi * mp.nint(mp.im(prev - v) / (2 * mp.pi))
        prev = v
    return prev


def Qf(w):
    _, _, _, dl = B.sigma(w, True)
    return dl * (w - L) * (w - Lb) / (L - Lb)


def Lamf(w):
    sw, _, _ = B.sigma(w)
    return sw / (zeta(w) * (L - Lb))


mp.mp.dps = 30
probe = []
worst = 0
for w in [mp.mpc(2, 2), mp.mpc(5, "1.6"), mp.mpc(20, 3), mp.mpc(1, "0.3"), mp.mpc(-3, "0.5"), mp.mpc("2.7", "0.05"),
          mp.mpc(60, "1.6"), mp.mpc(0, 3), mp.mpc("1.0", "0.01")]:
    d = dhyp(w)
    lq = stetig_log(Qf, w)
    lam = stetig_log(Lamf, w)
    g1 = abs(lq) / d
    g2 = abs(lam + mp.log(1 - abs(zeta(w)) ** 2)) / d
    worst = max(worst, g1, g2)
    probe.append([mp.nstr(w, 5), mp.nstr(d, 5), mp.nstr(g1, 4), mp.nstr(g2, 4)])
    print(f"w = {mp.nstr(w, 5):>14}: d = {mp.nstr(d, 5):>8}  |log Q|/d = {mp.nstr(g1, 4)}  |Λ+log(1−|ζ|²)|/d = {mp.nstr(g2, 4)}")
print("grösstes Verhältnis (muss ≤ 1 sein):", mp.nstr(worst, 4))
out["C"] = probe

# ---------------- Teil D: Konstanten auf D_r
mp.mp.dps = 30
# Bereich U3: Im u3 ∈ (π/2, π), Re u3 = ln|u2| ≥ X0 = ln ln(1/(A r (1+e1)))
X0 = mp.log(mp.log(1 / (A * r * (1 + e1))))
print("X0 = min Re u3 =", mp.nstr(X0, 5))


def W(z):
    l1 = mp.log(1 / (A * (1 - e1) * abs(z)))
    return mp.log(l1 + mp.pi) + mp.pi + abs(L)


def BLam(Wv):
    return mp.log(2 + 2 * Wv ** 2 / (mp.pi * ImL)) + mp.log(Wv ** 2 / (2 * mp.pi * ImL))


# sup von e^{±d} K(w) über U3 (Gitter, X ≥ X0)
H1 = mp.mpf(0)
for X in [X0 + mp.mpf(j) / 10 for j in range(0, 400)]:
    for Y in [mp.pi / 2 + (mp.pi / 2) * mp.mpf(i) / 20 for i in range(21)]:
        w = mp.mpc(X, Y)
        K = 2 * ImL / (abs(L) * abs(w - L) * abs(w - Lb))
        H1 = max(H1, mp.exp(dhyp(w)) * K)
print("H1 = sup_U3 e^d K =", mp.nstr(H1, 5), " (Grenzwert für Re→∞: 4/(π|L|) =", mp.nstr(4 / (mp.pi * abs(L)), 5), ")")
z1 = mp.mpc(0, "0.05")
th1 = T.theta(z1)
CL = 4 / abs(L)
# |θ(z)| ≤ |θ(z1)| + r + |z1| + (2·|Log ζ|max + BΛ(W(z)) + BΛ(W(z1)))/|L|;  |Log ζ| ≤ ln(1/(1−q)), q = sup|1−ζ|
qz = 2 * ImL / abs(mp.mpc(X0, mp.pi / 2) - Lb)
logz = -mp.log(1 - qz)
Cth = mp.mpf(0)
for e in range(1, 400):
    t = mp.mpf(10) ** (-mp.mpf(e) / 4)
    if t >= r:
        continue
    zt = mp.mpc(0, t)
    Wv = W(zt)
    bound = abs(th1) + r + abs(z1) + (2 * logz + BLam(Wv) + BLam(W(z1))) / abs(L)
    Cth = max(Cth, bound - CL * mp.log(Wv))
print(f"Majorante: |θ(z)| ≤ C_θ + (4/|L|) ln W(|z|),  C_θ = {mp.nstr(Cth, 5)}, 4/|L| = {mp.nstr(CL, 5)}")
# Ableitung: |1+θ'| ≤ (1+es)A·H1/(|u0||u1||u2|),  |u0| ≥ 1 − A r(1+e0), |u1| ≥ A(1−e1)|z|, |u2| ≥ ln(1/(A(1+e1)|z|))
Cplus = (1 + es) * A * H1 / ((1 - A * r * (1 + e0)) * A * (1 - e1))
print(f"|1+θ'(z)| ≤ C₊/(|z| ln(1/(A(1+ε1)|z|))),  C₊ = {mp.nstr(Cplus, 5)}")

# Vergleich mit dem wirklichen θ längs z = t e^{iφ}
k1 = int(k1)
vergl = []
mp.mp.dps = 60
for t in [mp.mpf("0.05"), mp.mpf("1e-3"), mp.mpf("1e-8"), mp.mpf("1e-20"), mp.mpf("1e-40")]:
    for phi in [mp.mpf("0.02"), mp.mpf("0.5"), mp.mpf("0.98")]:
        z = t * mp.expjpi(phi)
        s, sd = sexp0(z, True)
        u0 = s
        u1 = mp.log(u0)
        u2 = mp.log(u1)
        u3 = mp.log(u2)
        hv, hd = B.h_mit_ableitung(u3)
        th = 3 - z + hv + 2j * mp.pi * k1 / L
        w1 = sd * hd / (u0 * u1 * u2)
        Wv = W(z)
        Mth = Cth + CL * mp.log(Wv)
        Mw = Cplus / (t * mp.log(1 / (A * (1 + e1) * t)))
        vergl.append([mp.nstr(t, 3), mp.nstr(phi, 3), mp.nstr(abs(th), 6), mp.nstr(Mth, 5), mp.nstr(abs(w1) * t, 5),
                      mp.nstr(Mw * t, 5), mp.nstr(mp.im(u3), 4), mp.nstr(mp.arg(w1), 4)])
print(" t        φ/π   |θ|        Majorante  t|1+θ'|   t·Majorante  Im u3   arg(1+θ')")
for v in vergl:
    print("  ".join(f"{x:>9}" for x in v))
# Probe gegen die θ-Reihe an einem Punkt mit Im z ≥ 0.035
mp.mp.dps = 50
zt = mp.mpf("0.06") * mp.expjpi(mp.mpf("0.6"))
s = sexp0(zt)
u3 = mp.log(mp.log(mp.log(s)))
print("Probe θ-Reihe vs. 3 − z + h(u3) + 2πik/L bei z =", mp.nstr(zt, 5), ":",
      mp.nstr(abs(T.theta(zt) - (3 - zt + B.h(u3) + 2j * mp.pi * k1 / L)), 3))
out["D"] = {"X0": mp.nstr(X0, 6), "H1": mp.nstr(H1, 6), "C_theta": mp.nstr(Cth, 6), "4/|L|": mp.nstr(CL, 6),
            "C_plus": mp.nstr(Cplus, 6), "k1": k1, "vergleich": vergl}
(HIER / "daten").mkdir(exist_ok=True)
(HIER / "daten" / "l1_schranken.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")

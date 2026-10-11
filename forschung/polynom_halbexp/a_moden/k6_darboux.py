"""Strang A: lokale Entwicklung von f_K an L (aus Schröder φ, Ψ und Knesers κ_m) und Darboux-Abgleich mit den Taylor-Koeffizienten.

f_K(z) = φ(√L σ W(q)),  σ = Ψ(z), q = σ^γ = e^{2πiu}, W = e^{Lε(q)}  (k5_singularitaet.py)
       = f_reg(z) + Σ_{n≥1} S_n(t),  S_n(t) = t^{nγ} Q(t)^n F_n(t),  t = z − L,  Q = (σ/t)^γ,
       F_n(t) = Σ_k φ_k (√L σ)^k [q^n] e^{k L ε(q)}.
Schreibweise der Leitung: f_K − L = √L t (1 + Σ c_m t^m) + Σ_n d_n t^{1+nγ} (1 + Σ_j e_{n,j} t^j),
Zweig: t^α := (−L)^α (1 − z/L)^α (Hauptzweige), gültig in D nahe L.
Darboux: [z^k] (1 − z/L)^α = L^{−k} Γ(k − α) / (Γ(−α) Γ(k + 1));  a_k = 2 Re Σ_{n,j} d_n e_{n,j} (−L)^{α_{n,j}} [z^k](1 − z/L)^{α_{n,j}}.
"""
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402

HIER = Path(__file__).resolve().parent
GEM = HIER.parent / "gemeinsam"
DATEN = HIER / "daten"
mp.mp.dps = 70
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
SL = mp.exp(L / 2)
G = 2j * mp.pi / L
LOGML = mp.log(-L)
d = json.loads((DATEN / "singularitaet_theta.json").read_text(encoding="utf-8"))
cx = lambda p: mp.mpc(mp.mpf(p[0]), mp.mpf(p[1]))  # noqa: E731
phi = [cx(p) for p in d["phi"]]
psi = [cx(p) for p in d["psi"]]
kap = {int(k): cx(v) for k, v in d["kappa"].items()}
NT = 26          # Reihenlänge in t
NQ = 4           # Ordnung in q


def mul(a, b, n):
    return [mp.fsum(a[i] * b[k - i] for i in range(k + 1)) for k in range(n)]


def exp_r(g, n):
    E = [mp.exp(g[0])] + [mp.mpc(0)] * (n - 1)
    for k in range(1, n):
        E[k] = mp.fsum(j * g[j] * E[k - j] for j in range(1, k + 1)) / k
    return E


def log1p_r(a, n):
    """log(a(t)) für a_0 = 1."""
    lg = [mp.mpc(0)] * n
    for k in range(1, n):
        lg[k] = a[k] - mp.fsum(j * lg[j] * a[k - j] for j in range(1, k)) / k
    return lg


# ε(q) aus κ (gleiche Fixpunktiteration wie k5)
eps = [mp.mpc(0)] * (NQ + 1)
for _ in range(NQ + 3):
    neu = [mp.mpc(0)] * (NQ + 1)
    for m in range(1, NQ + 1):
        ex = exp_r([2j * mp.pi * m * x for x in eps], NQ + 1)
        fak = [(1 if k == 0 else 0) - (-1) ** m * ex[k] for k in range(NQ + 1)]
        for k in range(NQ + 1 - m):
            neu[k + m] += kap[m] * fak[k]
    eps = neu
# P_n(k) = [q^n] e^{kLε(q)} für k = 0..NT
Pn = {k: exp_r([k * L * x for x in eps], NQ + 1) for k in range(NT + 1)}
sig = psi[:NT] + [mp.mpc(0)] * max(0, NT - len(psi))       # σ(t) = Ψ(L + t)
Qlog = log1p_r(sig[1:] + [mp.mpc(0)], NT)                  # log(σ/t)
Q = exp_r([G * x for x in Qlog], NT)                        # (σ/t)^γ
F = {n: [mp.mpc(0)] * NT for n in range(NQ + 1)}
pot = [mp.mpc(1)] + [mp.mpc(0)] * (NT - 1)
for k in range(NT):
    for n in range(NQ + 1):
        c = phi[k] * SL ** k * Pn[k][n]
        for i in range(NT):
            F[n][i] += c * pot[i]
    pot = mul(pot, sig, NT)
# F_0 = f_reg (Kontrolle), S_n = Q^n F_n (ohne t^{nγ})
print(f"F_0: f_reg(L) − L = {mp.nstr(F[0][0] - L, 3)}, f_reg'(L) − √L = {mp.nstr(F[0][1] - SL, 3)}")
S = {}
for n in range(1, NQ + 1):
    Qn = [mp.mpc(1)] + [mp.mpc(0)] * (NT - 1)
    for _ in range(n):
        Qn = mul(Qn, Q, NT)
    S[n] = mul(Qn, F[n], NT)
    print(f"S_{n}: [t^0] = {mp.nstr(S[n][0], 3)} (soll 0), d_{n} = [t^1] = {mp.nstr(S[n][1], 25)}")
dn = {n: S[n][1] for n in S}
print(f"Kontrolle d_1 = 2 L^(3/2) κ_1 = {mp.nstr(2 * L * SL * kap[1], 25)}")
for n in (1, 2):
    print(f"e_{n},j (j = 1..5): " + ", ".join(mp.nstr(S[n][1 + j] / S[n][1], 10) for j in range(1, 6))
          + f";  |e_{n},j|^(1/j) bei j = 10, 20: {[mp.nstr(abs(S[n][1 + j] / S[n][1]) ** (mp.mpf(1) / j), 4) for j in (10, 20)]}")


# Darboux
def g(k, al):
    return mp.exp(mp.loggamma(k - al) - mp.loggamma(-al) - mp.loggamma(k + 1)) * L ** (-k)


def darboux(k, J=20, nmax=3):
    s = mp.mpc(0)
    for n in range(1, nmax + 1):
        for j in range(0, min(J + 1, NT - 1)):
            al = 1 + n * G + j
            s += S[n][1 + j] * mp.exp(al * LOGML) * g(k, al)
    return 2 * mp.re(s)


a = [mp.mpf(x) for x in json.loads((GEM / "kneser_halbexp_koeff_1000.json").read_text(encoding="utf-8"))["a"]]
fe = [mp.mpf(x) for x in json.loads((GEM / "kneser_halbexp_koeff_1000.json").read_text(encoding="utf-8"))["fehler_abs"]]
print("\nDarboux-Abgleich: a_k gegen 2Re Σ_{n≤3, j≤20} d_n e_{n,j} (−L)^α [z^k](1−z/L)^α")
print("  k   a_k                         rel.Abw.(n≤3,J=20)  rel.Abw.(nur n=1,j=0)  rel.Abw.(n=1, J=20)  Datenfehler rel.")
for k in (20, 40, 60, 100, 150, 200, 300, 400, 500, 600, 700, 800, 1000):
    p = darboux(k)
    p0 = darboux(k, J=0, nmax=1)
    p1 = darboux(k, J=20, nmax=1)
    print(f"{k:4d}  {mp.nstr(a[k], 18):>26}  {mp.nstr(abs(p - a[k]) / abs(a[k]), 3):>12}  {mp.nstr(abs(p0 - a[k]) / abs(a[k]), 3):>12}"
          f"  {mp.nstr(abs(p1 - a[k]) / abs(a[k]), 3):>12}  {mp.nstr(fe[k] / abs(a[k]), 3):>10}")
# Hüllkurve: a_k ≈ 2|D|/|Γ(−1−γ)| |L|^{−k} k^{−2−Re γ} cos(…)
D1 = dn[1] * mp.exp((1 + G) * LOGML)
amp = 2 * abs(D1 / mp.gamma(-1 - G))
print(f"\nLeitterm in (1 − z/L)-Form: D_1 = d_1 (−L)^(1+γ) = {mp.nstr(D1, 20)};  "
      f"a_k ≈ 2Re[D_1 k^(−2−γ) L^(−k) / Γ(−1−γ)],  2|D_1/Γ(−1−γ)| = {mp.nstr(amp, 12)}")
aus = {"beschreibung": ("Lokale Entwicklung von f_K = T_e(slog_e + 1/2) am Fixpunkt L (Strang A, k5/k6). t = z − L, t^α := (−L)^α(1 − z/L)^α "
                        "(Hauptzweige), γ = 2πi/L. f_K − L = √L t (1 + Σ_{m≥1} c_m t^m) + Σ_{n≥1} d_n t^{1+nγ}(1 + Σ_{j≥1} e_{n,j} t^j). "
                        "Analytischer Teil = reguläre (Schröder-)Halbiterierte f_reg = φ(√L Ψ). κ_m = Fourier-Koeffizienten von Knesers "
                        "ϑ(u) = A(z) − u, u = log Ψ(z)/L (gleicher Zweig). d_1 = 2 L^{3/2} κ_1. Darboux: a_k = 2Re Σ d_n e_{n,j} (−L)^α "
                        "L^{−k}Γ(k−α)/(Γ(−α)k!) mit α = 1+nγ+j; geprüft gegen kneser_halbexp_koeff_1000.json."),
       "L": [mp.nstr(mp.re(L), 40), mp.nstr(mp.im(L), 40)], "gamma": [mp.nstr(mp.re(G), 40), mp.nstr(mp.im(G), 40)],
       "sqrtL": [mp.nstr(mp.re(SL), 40), mp.nstr(mp.im(SL), 40)],
       "c_m": [[mp.nstr(mp.re(F[0][m + 1] / SL), 30), mp.nstr(mp.im(F[0][m + 1] / SL), 30)] for m in range(1, NT - 1)],
       "kappa": {str(k): [mp.nstr(mp.re(v), 35), mp.nstr(mp.im(v), 35)] for k, v in kap.items() if k >= 0},
       "d_n": {str(n): [mp.nstr(mp.re(dn[n]), 30), mp.nstr(mp.im(dn[n]), 30)] for n in dn},
       "e_nj": {str(n): [[mp.nstr(mp.re(S[n][1 + j] / S[n][1]), 25), mp.nstr(mp.im(S[n][1 + j] / S[n][1]), 25)] for j in range(1, NT - 1)]
                for n in (1, 2, 3)},
       "D1_form_1_minus_z_durch_L": [mp.nstr(mp.re(D1), 30), mp.nstr(mp.im(D1), 30)]}
(GEM / "kneser_singularitaet_L.json").write_text(json.dumps(aus, indent=1, ensure_ascii=False), encoding="utf-8")
print("geschrieben: gemeinsam/kneser_singularitaet_L.json")

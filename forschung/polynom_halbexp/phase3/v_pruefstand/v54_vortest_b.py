"""Strang V — Vortest (b) (K BRETT 11:00): duale, gewichtete Parametrix der Eck-Pol-Antwort u_n = J_n⁻¹T_nk_L (k_L = 1/(L − x)), kubisch.

(1) ũ_n = T_n[s] + T_n[g_n], s ≈ Koeffizienten ≤ n von u_{Mn} (M = 4, bei n = 800: 3), g_n(x) = (C/d_n)V((x − L)/d_n) − c_s/(x − L) = (C/d_n)(V(ζ) − A/ζ), A = c_s/C
    (V(ζ) = ∫₀¹F_P(λ̃_Psζ)ds über die Integraldarstellung F_β(u) = Γ(β)⁻¹∫t^{β−1}e^{−t−ut^{−iα}}dt; im Potenzsektor ist V − A/ζ exponentiell klein).
    Koeffizienten von g_n per Cauchy/FFT auf |x| = ρ(1 − δ), δ = 0.01, 8n Punkte (nicht über die Taylor-Reihe von V). Gegenproben: ũ = 0, ũ = T_n[s].
(2) r̃_n = b − J̃_nũ_n (ρ-skaliert, b_k = ρ^k/L^{k+1}), r̃_n(t) = Σ r̃_k(t/ρ)^k.
(3) w_n(t) = |Λ_n(t)/Λ_n(L)|, ℓ_n = J̃_n^{−T}e_n; R_n(t) = |r̃_n(t)|·w_n(t) auf 𝒞 (Σ₁ 80°, Σ₂ −90° bis 0.05, Verbindungsbögen, Linse 0.03, Γout-Bogen).
(4) Lemma ρ: ρ_{m,n} = (L⁻¹r̃_n)_m/(L⁻¹b)_m (LU ohne Pivot von J̃_n), m ≤ n; dazu ρ_{m,n−1} mit r̃_{n−1}.
Aufruf: python v54_vortest_b.py n [Mfaktor]
"""
import cmath
import json
import math
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np

HIER = Path(__file__).resolve().parent
NN = int(sys.argv[1])
MF = int(sys.argv[2]) if len(sys.argv) > 2 else 4
sys.argv = ["v48", "kubisch", str(MF * NN), "300"]
sys.path.insert(0, str(HIER))
import v48_RL as V48  # noqa: E402
from flint import arb, acb, arb_mat, arb_poly, ctx  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
t0 = time.time()
mp.mp.dps = 30
Lm = mp.mpc(V48.L_RE, V48.L_IM)
rho_m = abs(Lm)
al = Lm / (4 * mp.pi)
a_m = 1 + 1j * al
mu = mp.sqrt(Lm)
bp = 0.5 + 1j * al
lamP = mp.mpc("0.603424", "1.165020")
lt = -lamP * rho_m ** 2 / Lm ** 2
Cc = -lt * mp.gamma(0.5 + 1j * al) / ((mu + 1 / mu) * mp.gamma(0.5 + 2j * al))
cs = -1 / (mu + 1 / mu)
Lc = complex(Lm)
rho = abs(Lc)
th = cmath.phase(Lc)
cP = [mp.gamma(bp - 1j * al * k) / (mp.factorial(k) * mp.gamma(bp)) for k in range(200)]


def V_ser(z):
    s, p, w = mp.mpc(0), mp.mpc(1), -lt * z
    for k, c in enumerate(cP):
        s += c * p / (k + 1)
        p *= w
    return s


def V_int(z):
    def f(t):
        w = lt * z * mp.power(t, -1j * al)
        g = (1 - mp.exp(-w)) / w if abs(w) > 1e-12 else 1 - w / 2
        return mp.power(t, bp - 1) * mp.exp(-t) * g
    return mp.quad(f, [0, 1, 10, 100, 1000, mp.inf]) / mp.gamma(bp)


# Asymptotik im Potenzsektor (Mellin–Barnes, Pole rechts): V = A/ζ − Σ_m Γ(s'_m)(−1)^m (λ̃ζ)^{−s'_m}/(m!·iα·(1 − s'_m)Γ(β)), s'_m = −(β + m)/(iα)
A_as = mp.gamma(bp + 1j * al) / (mp.gamma(bp) * lt)
s_pole = [-(bp + m) / (1j * al) for m in range(6)]


def V_as(z):
    w = lt * z
    v = A_as / z
    for m, sm in enumerate(s_pole):
        v -= mp.gamma(sm) * (-1) ** m * mp.exp(-sm * mp.log(w)) / (mp.factorial(m) * 1j * al * (1 - sm) * mp.gamma(bp))
    return v


def V(z):
    if abs(lt * z) < 4:
        return V_ser(z)
    deg = float(mp.degrees(mp.arg(z))) % 360
    # Potenzsektor (Kontrolle gegen das Integral: rel. Fehler ≤ 1e-11 für arg ζ ∈ [200°, 345°] ab |ζ| = 30, ≤ 1e-8 bei 195°/345°)
    if (200 <= deg <= 350 and abs(z) >= 25) or (190 <= deg < 200 and abs(z) >= 80):
        return V_as(z)
    return V_int(z)


# A aus der Anpassung: (C/d_n)·A/ζ = c_s/(x − L) ⇒ A = c_s/C; Kontrolle gegen Γ(β+iα)/(Γ(β)λ̃)
A_pol = cs / Cc
print(f"  A = c_s/C = {mp.nstr(A_pol, 10)}, Mellin-Residuum Γ(β+iα)/(Γ(β)λ̃) = {mp.nstr(A_as, 10)}", flush=True)
for zz in (26 * mp.expjpi(1.0), 26 * mp.expjpi(-0.5), 40 * mp.expjpi(-0.8), 26 * mp.expjpi(-0.1), 60 * mp.expjpi(0.97)):
    vi, va = V_int(zz), V_as(zz)
    print(f"  Kontrolle bei ζ = {mp.nstr(zz, 3)}: Integral {mp.nstr(vi, 12)}, Asymptotik {mp.nstr(va, 12)}, |Δ| = {mp.nstr(abs(vi - va), 3)}", flush=True)


def g_koeff(n, delta=0.01, Zmax=1e9):
    """ρ-skalierte Taylor-Koeffizienten (k ≤ n) von g_n = (C/d_n)(V(ζ) − A/ζ), FFT auf |x| = ρ(1 − δ)."""
    dn = rho_m * mp.power(n, -a_m)
    M = 8 * n
    vals = np.zeros(M, dtype=complex)
    nint = 0
    for j in range(M):
        x = rho_m * (1 - delta) * mp.expjpi(2.0 * j / M)
        z = (x - Lm) / dn
        if abs(z) > Zmax:
            continue                                     # Potenzsektor: V − A/ζ exponentiell klein (oben kontrolliert)
        vals[j] = complex(Cc / dn * (V(z) - A_pol / z))
        nint += 1
    c = np.fft.fft(vals) / M
    k = np.arange(n + 1)
    return c[:n + 1] * (1 - delta) ** (-k.astype(float)), nint


def J_mat(N, P):
    ctx.prec = P
    V48.N = N
    fk, mk = V48.modell(N)
    f = arb_poly(fk)
    rh = (arb(V48.L_RE) ** 2 + arb(V48.L_IM) ** 2).sqrt()
    rp = [rh ** k for k in range(N + 1)]
    A = arb_mat(N + 1, N + 1)
    col = arb_poly([1])
    for j in range(N + 1):
        for k in range(N + 1):
            ck = col[k] if k <= col.degree() else arb(0)
            mm = mk[k - j] if k >= j else arb(0)
            A[k, j] = ((ck + mm) * rp[k] / rp[j]).mid()
        col = (col * f).truncate(N + 1)
    return A, rh


def bvec(N, rh):
    La = acb(arb(V48.L_RE), arb(V48.L_IM))
    q, p = acb(rh) / La, acb(1) / La
    out = []
    for k in range(N + 1):
        out.append(p)
        p = p * q
    return out


# ---- (1) u_{Mn}
NM = MF * NN
P = int(0.5 * NM) + 300
JM, rh = J_mat(NM, P)
bM = bvec(NM, rh)
B2 = arb_mat(NM + 1, 2)
for k in range(NM + 1):
    B2[k, 0], B2[k, 1] = bM[k].real.mid(), bM[k].imag.mid()
X = JM.solve(B2, algorithm="approx")
s_all = np.array([complex(float(X[k, 0].mid()), float(X[k, 1].mid())) for k in range(NM + 1)])
print(f"[n = {NN}] u_{NM} gelöst ({P} Bit, {time.time() - t0:.0f}s)", flush=True)
Pn = int(0.5 * NN) + 300
Jn, rh = J_mat(NN, Pn)
bn = bvec(NN, rh)
ctx.prec = Pn


def residuum(Y, n):
    """r̃ = b − J̃_n Y (Y komplex, Länge n+1), als komplexe numpy-Liste (aus arb)."""
    Yr = arb_mat(n + 1, 2)
    for k in range(n + 1):
        Yr[k, 0], Yr[k, 1] = arb(Y[k].real), arb(Y[k].imag)
    if n == NN:
        Jsub = Jn
    else:
        Jsub = arb_mat(n + 1, n + 1)
        for i in range(n + 1):
            for j in range(n + 1):
                Jsub[i, j] = Jn[i, j]
    Z = Jsub * Yr
    return [acb(bn[k].real - Z[k, 0], bn[k].imag - Z[k, 1]) for k in range(n + 1)]


gn, nint = g_koeff(NN)
gm, _ = g_koeff(NN - 1)
print(f"  g_n-Koeffizienten (FFT, {nint} Punkte mit |ζ| ≤ 60 ausgewertet) ({time.time() - t0:.0f}s)", flush=True)
varianten = {
    "null": (np.zeros(NN + 1, dtype=complex), np.zeros(NN, dtype=complex)),
    "bulk": (s_all[:NN + 1], s_all[:NN]),
    "komposit": (s_all[:NN + 1] + gn, s_all[:NN] + gm),
}
# ---- (3) adjungierter Kopf n
e = arb_mat(NN + 1, 1)
e[NN, 0] = arb(1)
ell = Jn.transpose().solve(e, algorithm="approx")
ellc = [ell[k, 0] for k in range(NN + 1)]


def Lam(t):
    q = acb(rh) / acb(t.real, t.imag)
    s = acb(0)
    for c in reversed(ellc):
        s = s * q + c
    return s


LamL = Lam(Lc)
# Konturpunkte (obere Hälfte; J̃ reell, also r̃(t̄) = conj(r̃(t)) nur wenn ũ reell — ũ ist komplex ⇒ beide Hälften)
TS = math.radians(1.623)


def spirale(uA, r, Lx):
    sgn = 1 if Lx.imag > 0 else -1
    return Lx + r * cmath.exp(1j * sgn * (math.radians(uA + 71.46) + TS * math.log(r / rho)))


kont = []
for Lx in (Lc, Lc.conjugate()):
    tag = "L" if Lx.imag > 0 else "L̄"
    for r in np.geomspace(1e-6, 0.05, 25):
        kont.append((spirale(80, r, Lx), f"Σ1@{tag}"))
        kont.append((spirale(-90, r, Lx), f"Σ2@{tag}"))
for phi in np.linspace(-th + 0.04, th - 0.04, 120):
    kont.append((rho * cmath.exp(1j * phi), "Γout"))
for phi in np.linspace(th + 0.04, 2 * math.pi - th - 0.04, 160):
    p = rho * 1.03 * cmath.exp(1j * phi)
    kont.append((p, "Steg" if abs(math.degrees(phi) - 180) <= 30 else "Linse"))
wts = {}
for p, k in kont:
    wts[(p, k)] = float(abs(Lam(p) / LamL).mid())
print(f"  Gewichte w_n auf 𝒞 ({len(kont)} Punkte) ({time.time() - t0:.0f}s)", flush=True)
erg = {}
res = {}
for name, (Yn, Ym) in varianten.items():
    r = residuum(Yn, NN)
    rm = residuum(Ym, NN - 1) + [acb(0)]
    res[name] = (r, rm)
    rows = {}
    for p, k in kont:
        x = acb(p.real, p.imag) / acb(rh)
        sv = acb(0)
        for c in reversed(r):
            sv = sv * x + c
        val = float(abs(sv).mid())
        R = val * wts[(p, k)]
        if k not in rows or R > rows[k][0]:
            rows[k] = (R, val, wts[(p, k)])
    erg[name] = {k: v[0] for k, v in rows.items()}
    print(f"  ũ = {name}: max R_n (|r̃|, w) — " + ", ".join(f"{k}: {v[0]:.2e} ({v[1]:.1e}, {v[2]:.1e})" for k, v in sorted(rows.items())), flush=True)
# ---- (4) Lemma ρ: EINE Vorwärts-Elimination mit Spalten b und r̃_n, r̃_{n−1} aller Varianten
ctx.prec = Pn
namen = list(varianten)
ncol = 2 + 4 * len(namen)
Aw = arb_mat(Jn)
Bw = arb_mat(NN + 1, ncol)
for k in range(NN + 1):
    Bw[k, 0], Bw[k, 1] = bn[k].real.mid(), bn[k].imag.mid()
    for i, name in enumerate(namen):
        r, rm = res[name]
        Bw[k, 2 + 4 * i], Bw[k, 3 + 4 * i] = r[k].real.mid(), r[k].imag.mid()
        Bw[k, 4 + 4 * i], Bw[k, 5 + 4 * i] = rm[k].real.mid(), rm[k].imag.mid()
Y = []
for k in range(NN + 1):
    p_ = Aw[k, k]
    Y.append([complex(float(Bw[k, 2 * j].mid()), float(Bw[k, 2 * j + 1].mid())) for j in range(ncol // 2)])
    if k == NN:
        break
    m = NN + 1
    cvec = arb_mat(m, 1)
    rvec = arb_mat(1, m)
    for ii in range(k + 1, m):
        cvec[ii, 0] = (Aw[ii, k] / p_).mid()
    for jj in range(k + 1, m):
        rvec[0, jj] = Aw[k, jj]
    brow = arb_mat(1, ncol)
    for jj in range(ncol):
        brow[0, jj] = Bw[k, jj]
    Aw = (Aw - cvec * rvec).mid()
    Bw = (Bw - cvec * brow).mid()
Y = np.array(Y)
yb = Y[:, 0]
for i, name in enumerate(namen):
    rho_n = Y[:, 1 + 2 * i] / yb
    rho_m1 = Y[:NN, 2 + 2 * i] / yb[:NN]
    mx = int(np.argmax(np.abs(rho_n)))
    diff = np.abs(rho_n[:NN] - rho_m1)
    md = int(np.argmax(diff))
    erg[name]["rho_max"] = float(np.abs(rho_n).max())
    erg[name]["rho_max_m"] = mx
    erg[name]["drho_max_n"] = float(diff.max() * NN)
    prof = " ".join(f"m={m_}:{abs(rho_n[m_]):.1e}" for m_ in (NN // 8, NN // 4, NN // 2, 3 * NN // 4, NN - 10, NN))
    print(f"  {name}: Lemma ρ: max_m|ρ_mn| = {np.abs(rho_n).max():.3e} bei m = {mx}; |ρ_nn| = {abs(rho_n[NN]):.3e}; max_m|ρ_mn − ρ_m,n−1|·n = "
          f"{diff.max() * NN:.3e} bei m = {md}; Profil {prof}  ({time.time() - t0:.0f}s)", flush=True)
json.dump(erg, open(HIER / f"v54_vortest_b_n{NN}.json", "w", encoding="utf-8"), indent=0)

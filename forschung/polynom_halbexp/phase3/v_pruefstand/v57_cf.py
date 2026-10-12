"""Strang V — (Fam-koh)-Test (Auftrag Leitung 12:00; K §34.2/§34.4 (ii)): Carathéodory–Fejér-Repräsentant von p_N − f_K auf Int γ.

Nur der Keim zählt: jedes P̃ ≡ p_N mod x^{N+2} ist erlaubt. Mit k: 𝔻 → Int γ, k(0) = 0 (k(w) = ρ·w·e^{F(w)}, Theodorsen wie v12) gilt x^{N+2}A(Int γ) = w^{N+2}A(𝔻), also
  min_{P̃} sup_{Int γ}|P̃ − f_K| = min_g ‖(p_N − f_K)∘k + w^{N+2}g‖_∞ = ‖T_{N+2}(c)‖₂  (Carathéodory–Fejér/Schur),
c_m = [w^m](p_N − f_K)∘k, m ≤ N+1, T die untere Dreiecks-Toeplitz-Matrix. Schranken: max|c_m| ≤ CF ≤ Σ|c_m|.
Inkohärenter Vergleich (K): c ~ (Kopfabweichung)·e^{NF₀}, F₀ = F(0) = harmonisches Mittel der Linsenhöhe; kohärent = CF polynomial in N.
Daten: p_N (A bis 1000, v33 für 1500/2000), f_K: P-V-zertifiziert bis 640, darüber A (4000, Darboux, rel. ~1e-21).
Aufruf: python v57_cf.py N1,N2,… [S] [tau1] [p]
"""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HIER = Path(__file__).resolve().parent
NS = [int(x) for x in sys.argv[1].split(",")]
S_ = sys.argv[2] if len(sys.argv) > 2 else "0.03"
T1_ = sys.argv[3] if len(sys.argv) > 3 else "12"
P_ = sys.argv[4] if len(sys.argv) > 4 else "1.3"
sys.argv = ["v12", "kubisch", "10", S_, T1_, P_, "K"]
sys.path.insert(0, str(HIER))
import v12_konform_S as V12  # noqa: E402
from flint import acb, arb, acb_poly, ctx  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
t0 = time.time()
F, it, d = V12.theodorsen()
res, kmin, kmax, mono, dL = V12.pruefe_abbildung(F)
F0 = F[0].real
print(f"[Linse K-Form S = {S_}, τ₁ = {T1_}°, p = {P_}] Theodorsen {it} It., Randfehler {res:.1e}, |k|/ρ ∈ [{kmin:.4f}, {kmax:.4f}], F₀ = {F0:.5f} ({time.time() - t0:.0f}s)", flush=True)
GEM = HIER.parent.parent / "gemeinsam"
zert = json.loads((GEM / "kneser_halbexp_koeff_zert_640.json").read_text(encoding="utf-8"))["a"]
a4000 = json.loads((GEM / "kneser_halbexp_koeff_4000.json").read_text(encoding="utf-8"))["a"]
pA = json.loads((GEM / "pN_N100_bis_1000.json").read_text(encoding="utf-8"))["p"]


def fK(j):
    x = zert[j] if j <= 640 else a4000[j]
    return arb(x[0] if isinstance(x, list) else x)


def pN(N):
    if N <= 1000:
        return [arb(s) for s in pA[str(N)]]
    return [arb(s) for s in json.loads((HIER / "sg" / f"v33_pN_{N}.json").read_text(encoding="utf-8"))["p"]]


erg = {}
for N in NS:
    t1 = time.time()
    ctx.prec = 500
    try:
        p = pN(N)
    except (KeyError, FileNotFoundError):
        print(f"N = {N}: keine Daten")
        continue
    M = N + 2
    rho = abs(acb(arb("0.31813150520476413531265425158766451720351761387139986692237"), arb("1.33723570143068940890116214319371061253950213846051241687")))
    # δ̃_j = (p_j − a_j)ρ^j, j ≤ N+1
    dt = []
    rj = arb(1)
    for j in range(M):
        pj = p[j] if j <= N else arb(0)
        dt.append((pj - fK(j)) * rj)
        rj = rj * rho
    # E = exp(F) bis Grad N+1 (Rekursion E_m = (1/m)Σ k F_k E_{m−k}), q = w·E
    Fa = [acb(arb(F[k].real), arb(F[k].imag)) for k in range(M)]
    E = [Fa[0].exp()]
    for m in range(1, M):
        sacc = acb(0)
        for k in range(1, m + 1):
            sacc += k * Fa[k] * E[m - k]
        E.append(sacc / m)
    q = acb_poly([acb(0)] + E[:M - 1])
    c = acb_poly([dt[0]])
    qj = acb_poly([1])
    for j in range(1, M):
        qj = (qj * q)
        qj = acb_poly([qj[i] for i in range(min(M, qj.degree() + 1))])
        c = c + qj * dt[j]
    cm = np.array([complex(float(c[i].real.mid()), float(c[i].imag.mid())) if i <= c.degree() else 0j for i in range(M)])
    am = np.abs(cm)
    # CF = ‖T(c)‖₂ (untere Dreiecks-Toeplitz)
    from scipy.linalg import toeplitz
    T = toeplitz(cm, np.zeros(M, dtype=complex))
    cf = np.linalg.norm(T, 2)
    ink = abs(float((dt[N].mid()))) * math.exp(N * F0)
    mmax = int(np.argmax(am))
    erg[N] = dict(CF=cf, max_c=float(am.max()), sum_c=float(am.sum()), m_max=mmax, inkoh=ink, F0=F0)
    print(f"N = {N} ({time.time() - t1:.0f}s): CF = min sup_{{Int γ}}|P̃ − f_K| = {cf:.3e} (Schranken {am.max():.3e} … {am.sum():.3e}; max|c_m| bei m = {mmax}); "
          f"|δ̃_N|·e^(N F₀) = {ink:.3e}; |δ̃_N| = {abs(float(dt[N].mid())):.2e}, e^(N F₀) = {math.exp(N * F0):.2e}; Profil |c_m|: " +
          " ".join(f"m={m_}:{am[m_]:.1e}" for m_ in (0, M // 4, M // 2, 3 * M // 4, M - 5, M - 2, M - 1)), flush=True)
json.dump({str(k): v for k, v in erg.items()}, open(HIER / f"v57_cf_S{S_}.json", "w", encoding="utf-8"), indent=0)

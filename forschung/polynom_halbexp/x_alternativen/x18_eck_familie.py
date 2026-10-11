"""X18 — Familien-Stabilität des ECKTEILS: die grössten Singulärwerte σ₁…σ₄ von A_P⁻¹, A_P = T_N J(P)|_{P_N}
(Tangente J(P)v = v∘P + (P'∘P)v, ℓ²_ρ-skaliert), für
  f      : P = f_K (A_P = J_N),
  pN     : P = p_N (echte Lösung, gemeinsam/pN_N100_bis_1000.json) — KOHÄRENTE Störung, sup_{D̄}|P − f| = R,
  zh     : P = f + η, η zufälliges reelles Polynom mit Koeffizienten nur in k ∈ [0.9N, N] (hochfrequent), sup_{|x|=ρ}|η| = R,
  zt     : wie zh, aber k ≤ N/2 (tieffrequent),
  zv     : wie zh, alle k ≤ N,
mit R := sup|p_N − f_K| (aus x17). Frage: ändert eine Störung der Grösse R die Eckmoden (σ₁ ~ e^{0.012N}) nur relativ um O(N·R)
(multiplikativ stabil) oder um σ₁·N·R (additiv, Eckverstärkung)?
Zufallsvarianten: m_F = m + T_N[η'∘f] + T_N[(m'/f')·η] (erste Ordnung in η, exakt bis O(η²); f''∘f = m'/f'), damit der Abschneidefehler
von T_N f in f'∘f NICHT mitgemessen wird; Kontrolle ft: P = T_N f mit m = T_N[(T_Nf)'∘(T_Nf)] (nur Abschneidefehler).
Methode: arb, J_N⁻¹ einmal, Blockpotenz für die grössten Eigenwerte von G = J⁻ᵀJ⁻¹, Ritz-Werte in mpmath.
Sekanten-Varianten (𝒥_P statt Tangente): pS (P = p_N), sh/st/sv (P = f + η wie zh/zt/zv).
Aufruf: python x18_eck_familie.py N [Varianten, Standard f,pN,ft,zh,zt,zv]
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
N = int(sys.argv[1])
VAR = sys.argv[2].split(",") if len(sys.argv) > 2 else ["f", "pN", "ft", "zh", "zt", "zv"]
ctx.prec = int(1.8 * N) + 400
mp.mp.dps = 60
A_str = json.loads((HIER.parent / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text())["a"]
PD = json.loads((HIER.parent / "gemeinsam" / "pN_N100_bis_1000.json").read_text())["p"][str(N)]
Lm = mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337))
rho_f = float(abs(Lm))
rho = arb(mp.nstr(abs(Lm), 55))
rk = [arb(1)]
for k in range(1, N + 1):
    rk.append(rk[-1] * rho)


def mullow(a, b, n):
    return (a * b).truncate(n)


def inv_series(a, n):
    g = arb_poly([1 / a[0]])
    k = 1
    while k < n:
        k = min(2 * k, n)
        g = mullow(g, arb_poly([2]) - mullow(a.truncate(k), g, k), k)
    return g


def sup_rand(coef):
    c = np.array([float(x) for x in coef])
    z = rho_f * np.exp(1j * np.linspace(0, 2 * np.pi, 4001)[:-1])
    return np.max(np.abs(np.polynomial.polynomial.polyval(z, c)))


F_K = [arb(s) for s in A_str[: N + 1]]
F_P = [arb(s) for s in PD] + [arb(0)] * (N + 1 - len(PD))
R = sup_rand([F_P[k] - F_K[k] for k in range(N + 1)])  # sup|T_N(p_N − f)| auf |x| = ρ (Schwanz von f zählt hier nicht: J_N sieht nur T_N F)


LRe = float(mp.re(Lm))


def form(kind):
    """Deterministische Störungen (sup_{|x|=ρ} = R): e0 = 1 (verschiebt den Fixpunkt), e1 = x, eL = q(x) := x² − 2Re L·x + ρ² (Null in L, L̄;
    ändert P'(L) um 2i·Im L), eL2 = q² (Null 2. Ordnung in L), eQ = q·x^K mit K = N/2 (Null in L, mittelfrequent)."""
    q = np.array([rho_f ** 2, -2 * LRe, 1.0])
    g = {"e0": np.array([1.0]), "e1": np.array([0.0, 1.0]), "eL": q, "eL2": np.convolve(q, q),
         "eQ": np.concatenate([np.zeros(N // 2), q]), "e2": np.array([0.0, 0.0, 1.0]), "e3": np.array([0.0, 0.0, 0.0, 1.0]),
         "e4": np.array([0.0, 0.0, 0.0, 0.0, 1.0])}[kind]
    g = np.concatenate([g, np.zeros(N + 1 - len(g))])
    g *= R / sup_rand(g)
    return [arb(float(x)) for x in g]


def stoerung(kind, seed=1):
    rng = np.random.default_rng(seed)
    lo, hi = {"zh": (int(0.9 * N), N), "zt": (0, N // 2), "zv": (0, N)}[kind]
    g = np.zeros(N + 1)
    g[lo:hi + 1] = rng.standard_normal(hi - lo + 1) / rho_f ** np.arange(lo, hi + 1)
    g *= R / sup_rand(g)
    return [arb(float(x)) for x in g]


FK2 = arb_poly([arb(s) for s in A_str[: N + 3]])
EX2 = arb_poly([1 / arb.fac_ui(k) for k in range(N + 3)])
M_EX = mullow(EX2, inv_series(FK2.derivative(), N + 2), N + 2)        # m = eˣ/f_K', bis x^{N+1}
MP_F = mullow(M_EX.derivative(), inv_series(FK2.derivative(), N + 1), N + 1)   # m'/f' = f''∘f


def J_N(F, eta=None, sekante=False):
    Fp = arb_poly(F)
    if F is F_K:
        m = M_EX.truncate(N + 1)
    elif sekante:
        # Sekante 𝒥_P: M_P = m + ½(f''∘f)(P − f) + O((P−f)²); mod x^{N+1} zählt nur P − T_N f (der Schwanz f_hi ist O(x^{N+1}))
        m = M_EX.truncate(N + 1) + mullow(MP_F, arb_poly(eta), N + 1) * arb(0.5)
    elif eta is not None:
        e = arb_poly(eta)
        de = e.derivative()
        comp = arb_poly([0])
        FN = arb_poly(F_K)
        for k in range(de.degree(), -1, -1):
            comp = mullow(comp, FN, N + 1) + arb_poly([de[k]])
        m = M_EX.truncate(N + 1) + comp + mullow(MP_F, e, N + 1)
    else:
        d = Fp.derivative()
        m = arb_poly([0])
        for k in range(d.degree(), -1, -1):
            m = mullow(m, Fp, N + 1) + arb_poly([d[k]])
    A = arb_mat(N + 1, N + 1)
    pw = arb_poly([1])
    for j in range(N + 1):
        for k in range(N + 1):
            A[k, j] = (pw[k] + (m[k - j] if k >= j else 0)) * rk[k] / rk[j]
        pw = mullow(pw, Fp, N + 1)
    return A


def top_sigma(A, b=5, it=40):
    Ai = A.inv()
    G = Ai.transpose() * Ai
    rng = np.random.default_rng(7)
    X = arb_mat([[arb(float(x)) for x in row] for row in rng.standard_normal((N + 1, b))])
    for _ in range(it):
        Y = G * X
        # Orthonormierung über Cholesky von YᵀY (klein, mpmath)
        S = Y.transpose() * Y
        Sm = mp.matrix([[mp.mpf(S[i, j].mid().str(50, radius=False)) for j in range(b)] for i in range(b)])
        try:
            Lc = mp.cholesky(Sm)
        except ValueError:
            return [mp.nan] * 4
        Li = Lc ** -1
        T = arb_mat([[arb(mp.nstr(Li[j, i], 50)) for j in range(b)] for i in range(b)])  # (L⁻¹)ᵀ
        X = Y * T
    H = X.transpose() * G * X
    Hm = mp.matrix([[mp.mpf(H[i, j].mid().str(50, radius=False)) for j in range(b)] for i in range(b)])
    ew = sorted(mp.eigsy(Hm)[0], reverse=True)
    return [mp.sqrt(e) for e in ew[:4]]


t0 = time.time()
print(f"N = {N}: R = sup|T_N(p_N − f_K)| auf |x|=ρ = {R:.3e}, N·R = {N * R:.3e}", flush=True)
erg = {}
for v in VAR:
    if v == "pS":
        A = J_N(F_P, [F_P[k] - F_K[k] for k in range(N + 1)], sekante=True)
    elif v[0] == "s" and v[1:] in ("e0", "e1", "eL", "eL2", "eQ", "e2", "e3", "e4"):
        eta = form(v[1:])
        A = J_N([F_K[k] + e for k, e in enumerate(eta)], eta, sekante=True)
    elif v[:2] in ("sh", "st", "sv"):
        eta = stoerung("z" + v[1], int(v[2:]) if len(v) > 2 else 1)
        A = J_N([F_K[k] + e for k, e in enumerate(eta)], eta, sekante=True)
    elif v in ("f", "pN", "ft"):
        F = F_K if v == "f" else F_P if v == "pN" else list(F_K)
        A = J_N(F)
    else:
        eta = stoerung(v)
        A = J_N([F_K[k] + e for k, e in enumerate(eta)], eta)
    s = top_sigma(A)
    erg[v] = s
    rel = ""
    if v != "f" and "f" in erg:
        rel = "  rel. Abw.: " + ", ".join(mp.nstr((a - b) / b, 3) for a, b in zip(s, erg["f"]))
        rel += f"  [geteilt durch R: " + ", ".join(mp.nstr((a - b) / b / R, 3) for a, b in zip(s, erg["f"])) + "]"
    print(f"  {v:3s}: σ₁…σ₄(A_P⁻¹) = " + ", ".join(mp.nstr(x, 14) for x in s) + rel + f"  ({time.time() - t0:.0f}s)", flush=True)
if len(VAR) > 1:
    (HIER / f"x18_N{N}.json").write_text(json.dumps({"N": N, "R": R, "sigma": {k: [str(x) for x in v] for k, v in erg.items()}}), encoding="utf-8")

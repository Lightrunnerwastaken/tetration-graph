"""Strang V — explizite Bilanz „(S_γ) mit Konstante B ⇒ (H₀) ⇒ Ballradius R_N ⇒ Familienbedingung X4a“ (Auftrag Leitung 23:20, Punkt 1).

Kette (Ks Satz K1 mit Lemma K1″(a), Satz P-I, Lemma X8′/Satz X4a, H1/H2):
  E_N(η,z) = Π_Nk_η(z) − k_η(z), z ∈ D̄, η auf Schlitzen S_L, S_L̄ und Niveaulinie Γ_d = {g_γ = d}.
  K1″(a): ‖Π_Nk_η‖_{H²(D_ρ)} ≤ B(1 + M_γ)/dist(η, K_γ);  sup_{D̄}|P| ≤ √(N+1)‖P‖_{H²(D_ρ)}.
  (I)  dist(η, K_γ) ≥ 1/N: |E_N| ≤ [√(N+1)·N·B(1+M_γ) + (N+1)/ρ]·e^{(N+2)Lip/N}·e^{−(N+2)g_γ(η)} + (ρ/|η|)^{N+1}/(|η| − ρ)
       (Bernstein–Walsh in Ĉ∖K_γ^{1/N}, g_{1/N} ≥ g_γ − Lip/N, y_η = O(η^{−N−2}));
  (II) dist < 1/N:  |E_N| ≤ √(N+1)B(1+M_γ)/dist + 1/(|η| − ρ).
  P-I: sup_{D̄}|Π_Nf − f| ≤ 2·(1/2π)∫_0^{u₀}|[f](L(1+u))|·sup|E_N|·ρ du + (1/2π)M_d∮_{Γ_d} sup|E_N| |dη|.
  Sprung |[f](L(1+u))| ≤ C_J u^{5.447}, C_J = 0.95 für u ≤ 0.2 [num., A k25 / P]; d so gewählt, dass u₀(d) ≤ 0.2.
  Für die Familie P ∈ 𝒮_N ∩ B_N gilt dieselbe Schranke mit C_R^P = 2C_R (X4a) — hier konservativ Faktor 2 auf den B-Anteil.
  R_N = K₀·sup|Π^P_Nf − f| (K₀ = 6, H2), θ_N = A_θ·N·R_N ≤ ½ (X4a, eckfreier Teil der Familie), A_θ aus Lemma X8′ mit c₀ = 0.641 ((a″)_∞, V).
Geometrie aus v30 (g_γ-Tabelle längs des Schlitzes, Lip, M_γ, c_d, Niveaulinien). M_d = max|f_K| auf Γ_d über die Funktionalgleichungen.
Ausgabe: R_N = K₀(B·A_N + C_N), θ_N, und B_max(N₀) = größtes B mit θ_N ≤ ½ für alle N ∈ [N₀, 5000].
Aufruf: python v31_N0_bilanz.py S d
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad

HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER / "m_kopie"))
import mk_kneser as MK  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
S = sys.argv[1]
d = sys.argv[2]
G = json.loads((HIER / f"v30_Kgamma_S{S}.json").read_text(encoding="utf-8"))
rho, th = MK.RHO, MK.TH
L = MK.L
lip, Mg, cd = G["lip_x"], G["M_gamma"], G["c_d"]
if len(sys.argv) > 3:
    Mg = float(sys.argv[3])          # rigorose Schranke aus v32 (2.40 statt 2.3376)
uu, gg = np.array(G["g_tabelle"]["u"]), np.array(G["g_tabelle"]["g"])
lv = G["niveaulinien"][d]
u0, laenge = lv["u0"], lv["laenge"]
CJ, SIG = 0.95, 5.447
K0 = 6.0
c0 = 0.641
minfp, max_m_fp, M2 = 0.3456, 2.1176, 0.8183
assert u0 <= 0.2 + 1e-9, f"u₀ = {u0} > 0.2: Sprungschranke nicht verfügbar"


def A_theta(N):
    return (5 / rho) / minfp + max_m_fp / (rho * c0) + M2 / (c0 * N)


def g_u(u):
    """g_γ(L(1+u)) aus der Tabelle (log-log-Interpolation); unterhalb der Tabelle konservativ 0."""
    if u < uu.min():
        return 0.0
    return float(np.exp(np.interp(math.log(u), np.log(uu), np.log(gg))))


# M_d: |f_K| auf Γ_d
pts = np.array([complex(a, b) for a, b in lv["punkte"]])
vals = []
for x in pts:
    if abs(x) <= rho * 0.999:
        v = MK.F_innen(np.array([x / rho]))[0] * rho
    elif abs(np.log(x)) < rho * 0.97:
        v = np.exp(MK.F_innen(np.array([np.log(x) / rho]))[0] * rho)
    elif x.real < rho * math.cos(th) * 0.97:
        v = np.log(MK.F_innen(np.array([np.exp(x) / rho]))[0] * rho)
    else:
        v = np.nan
    vals.append(abs(v))
vals = np.array(vals)
nicht = int(np.sum(np.isnan(vals)))
Md = float(np.nanmax(vals))
Md_eff = max(Md, 2.5)          # nicht auswertbare Punkte liegen nahe L (|f| ≈ |L| ≈ 1.4–1.6) — Schranke 2.5 dort angenommen
print(f"[S = {S}, d = {d}] Lip = {lip:.4f}, M_γ = {Mg:.4f}, c_d = {cd:.4f}, u₀ = {u0:.4f}, |Γ_d| = {laenge:.3f}, M_d = {Md:.3f} "
      f"({nicht} von {len(vals)} Punkten nahe L nicht ausgewertet, dort |f| ≤ 2.5 angenommen)", flush=True)


def terme(N):
    """(A_N, C_N): sup|Π_Nf − f| ≤ B·A_N + C_N (Familie: Faktor 2 auf den B-Anteil)."""
    sq = math.sqrt(N + 1)
    ePre = math.exp((N + 2) * lip / N)

    def EN_B(u):           # Koeffizient von B in sup|E_N| auf dem Schlitz
        dist = cd * rho * u
        if dist >= 1.0 / N:
            return sq * N * (1 + Mg) * ePre * math.exp(-(N + 2) * g_u(u))
        return sq * (1 + Mg) / dist

    def EN_C(u):
        dist = cd * rho * u
        eta_abs = rho * (1 + u)
        if dist >= 1.0 / N:
            return (N + 1) / rho * ePre * math.exp(-(N + 2) * g_u(u)) + (rho / eta_abs) ** (N + 1) / (eta_abs - rho)
        return 1.0 / (eta_abs - rho)

    pre = 2 * rho / (2 * math.pi)
    A_s = pre * quad(lambda u: CJ * u ** SIG * EN_B(u), 0, u0, limit=400, points=[1 / (N * cd * rho)])[0]
    C_s = pre * quad(lambda u: CJ * u ** SIG * EN_C(u), 0, u0, limit=400, points=[1 / (N * cd * rho)])[0]
    dd = float(d)
    rmin = rho + lv["min_abs_minus_rho"]
    A_c = Md_eff * laenge / (2 * math.pi) * sq * N * (1 + Mg) * ePre * math.exp(-(N + 2) * dd)
    C_c = Md_eff * laenge / (2 * math.pi) * ((N + 1) / rho * ePre * math.exp(-(N + 2) * dd) + (rho / rmin) ** (N + 1) / (rmin - rho))
    return 2 * (A_s + A_c), C_s + C_c, (A_s, A_c)


print(" N    | A_N (·B)   | C_N        | Schlitz/Kontur (B-Teil) | θ_N/B (≈)  | B_max(N)", flush=True)
bmax = {}
for N in [200, 250, 300, 318, 350, 400, 500, 700, 1000, 1500, 2000, 3000, 5000]:
    A, C, (As, Ac) = terme(N)
    At = A_theta(N)
    # θ_N = A_θ·N·K₀·(B·A + C) ≤ ½
    b = (0.5 / (At * N * K0) - C) / A
    bmax[N] = b
    print(f" {N:4d} | {A:.3e} | {C:.3e} | {As:.2e} / {Ac:.2e} | {At * N * K0 * A:.3e} | {b:.4g}", flush=True)
for N0 in (318, 400):
    b0 = min(v for k, v in bmax.items() if k >= N0)
    print(f"⇒ B_max(N₀ = {N0}) = {b0:.4g}  (gemessen ‖Ĵ⁻¹‖ = 1.2503 Ks Linse / 1.1073 S = 0.05)", flush=True)

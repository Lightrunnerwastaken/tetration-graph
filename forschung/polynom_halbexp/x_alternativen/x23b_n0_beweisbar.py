"""X23b — N₀-Bilanz mit BEWEISBAREN Konstanten (Meilenstein 6): Lemma X7 explizit (C_X7 = 6.28 ⇒ C_R = 12.7, κ_R = 1),
A_θ beweisbar (19.7 bei c₀ = 0.54, Zeilen bis 5N; 15.5 bei 3N), Konturvarianten:
  K1: Niveaulinie d = 0.3 (Konstante 17.1), braucht Sprungschranke am Schlitz bis u₀(0.3) = 0.451;
  K2: Niveaulinie d = 0.1 (37.9), Sprung nur bis u = 0.151 nötig (heute belegt bis 0.2);
  K3: Hybrid — Niveaulinie d = 0.3 fern von L, Schlitz nur bis u = 0.2 (g = 0.134), dort Eck-Rate des Schlitzes ε_S ≈ 0.89 (W1: g_{K_T}/g_W̄ = 0.89–0.97)
      ⇒ Zusatzterm (2M/π)·C_E N^{κ_E}e^{−ε_S·0.134·N}·log(u₀/0.2)/c_d mit dem groben Sprung |Δf| ≤ 2M_d.
θ_N = A_θ·N·R_N ≤ ½; R_N/K₀ = 2·17.7·(2C_R)·N^{κ_R−4.2877} + 2·12.7·ε^{−4.8094}·C_E·N^{κ_E−4.8094} + Kontur.
"""
import numpy as np

S_K, e_K, S_W, e_W, K0 = 17.7, 4.2877, 12.7, 4.8094, 6


def RN(N, C_R, k_R, C_E, k_E, eps, kont):
    rest = 2 * S_K * 2 * C_R * N ** (k_R - e_K)
    eck = 2 * S_W * eps ** (-e_W) * C_E * N ** (k_E - e_W)
    return K0 * (rest + eck + kont(N, C_E, k_E, eps))


def K1(N, C_E, k_E, eps):
    return 17.1 * C_E * N ** k_E * np.exp(-eps * 0.3 * N)


def K2(N, C_E, k_E, eps):
    return 37.9 * C_E * N ** k_E * np.exp(-eps * 0.1 * N)


def K3(N, C_E, k_E, eps):
    zusatz = (2 * 4.25 / np.pi) * C_E * N ** k_E * np.exp(-0.89 * 0.134 * N) * np.log(0.451 / 0.2) / 0.933
    return K1(N, C_E, k_E, eps) + zusatz


def N0(A, **kw):
    for N in range(20, 10 ** 6):
        if A * N * RN(N, **kw) <= 0.5:
            return N
    return None


basis = dict(C_R=12.7, k_R=1.0, C_E=1.0, k_E=2.0, eps=0.69)
print("Kontur × A_θ (C_R = 12.7, κ_R = 1; Ecke C_E = 1, κ_E = 2, ε = 0.69):")
for name, kont in (("K1 d=0.3", K1), ("K2 d=0.1", K2), ("K3 Hybrid", K3)):
    row = []
    for A in (4.05, 13.1, 15.5, 19.7, 22.2):
        row.append(f"A_θ={A}: {N0(A, kont=kont, **basis)}")
    print(f"  {name}: " + "; ".join(row))
print("Eck-Konstanten bei K1, A_θ = 19.7 (beweisbar) bzw. 15.5 (Zeilen bis 3N):")
for kE in (1.0, 1.5, 2.0, 2.5):
    row = []
    for CE in (0.01, 0.1, 1, 10):
        b = dict(basis, C_E=CE, k_E=kE)
        row.append(f"C_E={CE}: {N0(19.7, kont=K1, **b)}/{N0(15.5, kont=K1, **b)}")
    print(f"  κ_E = {kE}: " + "; ".join(row))
print("Empfindlichkeit dN₀/dA_θ (K1, C_E = 1, κ_E = 2):")
for A in (2, 4.05, 8, 13.1, 19.7, 30, 50):
    print(f"  A_θ = {A}: N₀ = {N0(A, kont=K1, **basis)}")

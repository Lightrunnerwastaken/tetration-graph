"""X23 — N₀-Bilanz für KONSTANTEN.md: R_N und die Familien-Bedingung θ_N ≤ 1/2 (Satz X4a) in Szenarien.

R_N/K₀ = C'_R·N^{κ_R − e'} + C'_E·N^{κ_E − e_W} + Kontur(d),
  C'_R = 2·S_{K'}·C_R^P (Rest/eckfreier Teil in K', volle Green-Rate, Schwelle e' = 5.447·α' = 4.2877, S_{K'} = 17.7),
  C'_E = 2·S_{W̄}·ε_E^{−e_W}·C_E (Eckteil in W̄, Schwelle e_W = 4.8094, S_{W̄} = 12.7 (P)),
  Kontur = (M_d|Γ_d|/(2πδ_d))·C_E·N^{κ_E}·e^{−ε_E N d}.
θ_N = A_θ·N·R_N, A_θ = K_B(1 + max|m|/c₀) + ½M2/(c₀N), K_B = 0.76, max|m| = 2.3376, M2 = 0.8183 (x21).
C_R^P = 2·C_R (für 𝒥_P gilt (a″) nur mit c₀/2 bei θ ≤ 1/2).
"""
import numpy as np

S_K, e_K = 17.7, 4.2877
S_W, e_W = 12.7, 4.8094
K0 = 6
KB, mmax, M2 = 0.76, 2.3376, 0.8183


def RN(N, s):
    rest = 2 * S_K * (2 * s["C_R"]) * N ** (s["k_R"] - e_K)
    eck = 2 * S_W * s["eps"] ** (-e_W) * s["C_E"] * N ** (s["k_E"] - e_W)
    kontur = s["kontur"] * s["C_E"] * N ** s["k_E"] * np.exp(-s["eps"] * N * s["d"])
    return K0 * (rest + eck + kontur), K0 * rest, K0 * eck, K0 * kontur


def theta(N, s):
    A = KB * (1 + mmax / s["c0"]) + 0.5 * M2 / (s["c0"] * N)
    return A * N * RN(N, s)[0]


szen = {
    "S1 heute bewiesen (P-II: |w| ≤ eN^{3/2}/(c_fρ) ⇒ κ_R = 2, C_R = 10.6; Ecke wie gemessen C_E = 1, κ_E = 2; d = 0.1)":
        dict(C_R=10.6, k_R=2.0, C_E=1.0, k_E=2.0, eps=0.69, c0=0.54, d=0.1, kontur=37.9),
    "S2 P-II scharf (|w| ≤ 1.17√M ⇒ κ_R = 1, C_R = 2.9); Ecke C_E = 1, κ_E = 2; d = 0.1":
        dict(C_R=2.9, k_R=1.0, C_E=1.0, k_E=2.0, eps=0.69, c0=0.54, d=0.1, kontur=37.9),
    "S3 wie S2, aber Kontur d = 0.3 (x24: 17.1)":
        dict(C_R=2.9, k_R=1.0, C_E=1.0, k_E=2.0, eps=0.69, c0=0.54, d=0.3, kontur=17.1),
    "S4 wie S3, Ecke schlechter: C_E = 1, κ_E = 3":
        dict(C_R=2.9, k_R=1.0, C_E=1.0, k_E=3.0, eps=0.69, c0=0.54, d=0.3, kontur=17.1),
    "S5 wie S3, c₀ = 0.4249 (Wert von M6 als Platzhalter)":
        dict(C_R=2.9, k_R=1.0, C_E=1.0, k_E=2.0, eps=0.69, c0=0.4249, d=0.3, kontur=17.1),
    "S6 gemessen (R̂(L) ≈ 0.33N^{1.1}, Rc(L) ≈ 0.01N^{1.6}, ε = 0.87, d = 0.3)":
        dict(C_R=0.33, k_R=1.1, C_E=0.01, k_E=1.6, eps=0.87, c0=0.73, d=0.3, kontur=17.1),
}
for name, s in szen.items():
    N0 = next((N for N in range(10, 200000) if theta(N, s) <= 0.5 and (s["d"] >= 0.3 or RN(N, s)[3] <= RN(N, s)[2] + RN(N, s)[1])), None)
    r193 = RN(193, s)
    print(f"{name}:\n   R_193 = {r193[0]:.2e} (Rest {r193[1]:.1e}, Ecke {r193[2]:.1e}, Kontur {r193[3]:.1e}); θ_193 = {theta(193, s):.2e}; "
          f"N₀(θ ≤ 1/2) = {N0}")
# Grenz-κ_E für N₀ ≤ 193 in S3
for kE in (1.5, 2.0, 2.2, 2.5, 3.0):
    for CE in (0.1, 1, 10):
        s = dict(szen["S3 wie S2, aber Kontur d = 0.3 (x24: 17.1)"], k_E=kE, C_E=CE)
        N0 = next((N for N in range(10, 10 ** 6) if theta(N, s) <= 0.5), None)
        print(f"   S3 mit κ_E = {kE}, C_E = {CE}: N₀ = {N0}")

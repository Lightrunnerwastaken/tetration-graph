"""Univalenz von g(q) = θ − a₀ = Σ_{n≥1} a_n qⁿ in der Kreisscheibe (Basis e, a_n aus etappe2, n ≤ 4096).

Gebraucht für die lineare Starrheit der Momentenfamilie (Bericht, Satz M):
  ist g auf |q| ≤ ρ univalent, so erzwingt die linearisierte Momentenfamilie p = konst. für jedes p,
  das im Streifen |Im z| < δ mit e^{−2πδ} < ρ holomorph ist.
Prüfungen (Gleitkomma, double):
  (1) hinreichendes Kriterium |g'(q)/a₁ − 1| < 1 aus Σ_{n≥2} n|a_n| rⁿ⁻¹ < |a₁|  →  Radius r_suff
  (2) Nullstellen von g' in |q| < ρ (Windungszahl von g' auf |q| = ρ)
  (3) Randkurve g(ρe^{it}) doppelpunktfrei (Schnitt aller Segmentpaare)
(2)+(3) ⇒ g univalent auf |q| ≤ ρ (Darboux). Abbruchfehler der Reihe: Σ_{n>4096}|a_n|ρⁿ.
"""
import json
from pathlib import Path

import numpy as np

HIER = Path(__file__).resolve().parent
d = json.loads((HIER.parent / "etappe2" / "daten" / "a_e_80_0.003_16384.json").read_text(encoding="utf-8"))
a = np.array([complex(float(x[0]), float(x[1])) for x in d["a"]])
n = np.arange(len(a))
a1 = abs(a[1])
print(f"{len(a)-1} Koeffizienten, |a1| = {a1:.6f}, |a2| = {abs(a[2]):.6f}, |a_4096| = {abs(a[-1]):.3e}")

# (1)
def summe(r):
    return np.sum(n[2:] * np.abs(a[2:]) * r ** (n[2:] - 1))
lo, hi = 0.0, 1.0
for _ in range(60):
    m = (lo + hi) / 2
    (lo, hi) = (m, hi) if summe(m) < a1 else (lo, m)
print(f"(1) hinreichendes Kriterium: g univalent für |q| < r_suff = {lo:.4f}  (δ_suff = ln(1/r)/2π = {np.log(1/lo)/(2*np.pi):.4f})")


def werte(rho, N):
    t = np.linspace(0, 2 * np.pi, N, endpoint=False)
    q = rho * np.exp(1j * t)
    # Horner
    g = np.zeros(N, complex); dg = np.zeros(N, complex)
    for k in range(len(a) - 1, 0, -1):
        g = g * q + a[k]
        dg = dg * q + k * a[k]
    return g * q, dg


def doppelpunkte(z):
    P = np.stack([z.real, z.imag], 1)
    A = P; Bp = np.roll(P, -1, 0)
    N = len(P)
    cnt = 0
    for i in range(N):
        j = np.arange(i + 2, N if i > 0 else N - 1)
        if len(j) == 0:
            continue
        p, r = A[i], Bp[i] - A[i]
        q, s = A[j], Bp[j] - A[j]
        den = r[0] * s[:, 1] - r[1] * s[:, 0]
        qp = q - p
        with np.errstate(divide="ignore", invalid="ignore"):
            tt = (qp[:, 0] * s[:, 1] - qp[:, 1] * s[:, 0]) / den
            uu = (qp[:, 0] * r[1] - qp[:, 1] * r[0]) / den
        cnt += int(np.sum((tt >= 0) & (tt <= 1) & (uu >= 0) & (uu <= 1) & (den != 0)))
    return cnt


print("(2)+(3)  rho | Rest Σ_{n>4096}|a_n|ρⁿ | Nullst. g' | Doppelpunkte Randkurve | min|g'|")
for rho in (0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995):
    N = 4000 if rho < 0.95 else 12000
    g, dg = werte(rho, N)
    wind = np.sum(np.diff(np.unwrap(np.angle(np.append(dg, dg[0]))))) / (2 * np.pi)
    rest = abs(a[-1]) * rho ** len(a) / (1 - rho)
    dp = doppelpunkte(g)
    print(f"   {rho:5.3f} | {rest:9.2e} | {wind:+.3f} | {dp} | {np.min(np.abs(dg)):.3e}", flush=True)

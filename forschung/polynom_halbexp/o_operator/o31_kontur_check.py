"""O31 — Probe der Konturverschiebung (o29): X(δ) = Re⟨Kδ, δ⟩_p für Polynome δ einmal direkt auf ∂D (Quadratur in φ),
einmal über die Mittellinien:  X = Re∫ W_b u(z_b+τ_b) conj(u(z_b)) dξ + Re∫ W_t u(z_t) conj(u(z_t+τ_t)) dξ,  u = δ(x(ζ))·J(ζ)^{1/2}.
Normierung: Werte auf ∂D mit dφ/2π ⇔ (1/(2πρ))·dξ-Integrale im Streifen. Ausserdem D(u) und ‖Bu‖² (Strenge-Kette X ≥ −D).
"""
import json
import os
import sys

import numpy as np

import o26_streifen as G

sys.stdout.reconfigure(encoding="utf-8")
L, Lb, th0, RHO = G.L, G.Lb, G.th0, abs(G.L)
d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "gemeinsam",
                                "kneser_halbexp_koeff_4000.json"), encoding="utf-8"))
a = np.array([float(x) for x in d["a"]])
kk = np.arange(len(a))
f = lambda x: np.polynomial.polynomial.polyval(x, a)
fp = lambda x: np.polynomial.polynomial.polyval(x, a[1:] * kk[1:])
thL = np.angle(L)


def xz(z):
    w = np.exp(z + 1j * th0)
    return (L - Lb * w) / (1 - w)


def Jz(z):
    # x = (L − L̄ w)/(1 − w), w = e^{z+iθ₀}: dx/dz = w (L − L̄)/(1 − w)²
    w = np.exp(z + 1j * th0)
    return w * (L - Lb) / (1 - w) ** 2


rng = np.random.default_rng(1)
ub = np.load("o29_unten.npz")
ob = np.load("o29_oben.npz")
xi = ub["xi"]
h = xi[1] - xi[0]


def sqrtJ(z, ref):
    s = np.sqrt(Jz(z))
    # stetiger Zweig relativ zu ref (gleiche Länge)
    flip = np.abs(s - ref) > np.abs(-s - ref)
    s[flip] = -s[flip]
    return s


for trial in range(4):
    N = 8
    c = (rng.standard_normal(N + 1) + 1j * rng.standard_normal(N + 1)) / (np.arange(N + 1) + 1) / RHO ** np.arange(N + 1)
    delta = lambda x: np.polynomial.polynomial.polyval(x, c)
    # direkt: X = ∫ p (δ∘f/m) conj(δ) dφ/2π, p = |m|², m = e^x/f'
    ph = np.linspace(0, 2 * np.pi, 20000, endpoint=False)
    x = RHO * np.exp(1j * ph)
    m = np.exp(x) / fp(x)
    Xd = np.mean(np.real(np.abs(m) ** 2 * delta(f(x)) / m * np.conj(delta(x))))
    # über die Kurven: u = δ(x(ζ)) J^{1/2}
    zb = xi + 1j * ub["eps"]
    zbt = zb + ub["tau"]
    ref_b = np.sqrt(Jz(xi + 0j))
    for i in range(1, len(ref_b)):
        if abs(ref_b[i] - ref_b[i - 1]) > abs(-ref_b[i] - ref_b[i - 1]):
            ref_b[i] = -ref_b[i]
    sb = sqrtJ(zb, ref_b)
    sbt = sqrtJ(zbt, ref_b)
    u_b = delta(xz(zb)) * sb
    u_bt = delta(xz(zbt)) * sbt
    Xb = np.sum(np.real(ub["W"] * u_bt * np.conj(u_b))) * h
    zt = xi + 1j * (np.pi + ob["eps"])
    ztt = zt + ob["tau"]
    ref_t = np.sqrt(Jz(xi + 1j * np.pi))
    for i in range(1, len(ref_t)):
        if abs(ref_t[i] - ref_t[i - 1]) > abs(-ref_t[i] - ref_t[i - 1]):
            ref_t[i] = -ref_t[i]
    # Zweig oben an unten anschliessen: entlang ξ = −12 vertikal stetig
    yy = np.linspace(0, np.pi, 200)
    col = np.sqrt(Jz(xi[0] + 1j * yy))
    for i in range(1, len(col)):
        if abs(col[i] - col[i - 1]) > abs(-col[i] - col[i - 1]):
            col[i] = -col[i]
    if col[0] * ref_b[0] == 0 or abs(col[0] - ref_b[0]) > abs(col[0] + ref_b[0]):
        col = -col
    if abs(ref_t[0] - col[-1]) > abs(ref_t[0] + col[-1]):
        ref_t = -ref_t
    st = sqrtJ(zt, ref_t)
    stt = sqrtJ(ztt, ref_t)
    u_t = delta(xz(zt)) * st
    u_tt = delta(xz(ztt)) * stt
    Xt = np.sum(np.real(ob["W"] * u_t * np.conj(u_tt))) * h
    Xs = (Xb + Xt) / (2 * np.pi * RHO)
    print(f"Probe {trial}: X direkt = {Xd:+.8f},  über Mittellinien = {Xs:+.8f}  (unten {Xb / (2 * np.pi * RHO):+.6f}, oben {Xt / (2 * np.pi * RHO):+.6f})")

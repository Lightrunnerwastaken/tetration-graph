"""O26 — Geometrie in der Streifenkoordinate ζ = log M(x) − iθ₀, M(x) = (x − L)/(x − L̄):
D = exakter Streifen {0 < Im ζ < π}, Γin = {Im ζ = 0}, Γout = {Im ζ = π}; L ↔ Re ζ = −∞, L̄ ↔ +∞; f ↔ g(ζ).
Längs Γin und Γout (parametrisiert durch ξ = Re ζ): Verschiebung ℓ_b(ξ) = g(ξ) − ξ, ℓ_t(ξ) = g(ξ+iπ) − (ξ+iπ), Gewichte
R = (g'/f')^{1/2}/m (H²-isometrische Normierung δ = u·J^{−1/2}, J = dx/dζ), p = |m|², Höhe der Linsenoberkante h(ξ) = π + Im ℓ_t.
Dazu das Symbol des „eingefrorenen“ Modells bei ξ₀ (Fourier u = ∫û e^{iσζ}):
  B(σ) = p_b + p_t|R_t|² e^{−2σ(π + Im ℓ_t)},  X(σ) = p_b Re(R_b e^{iσℓ_b}) + p_t e^{−σ(2π + Im ℓ_t)} Re(R̄_t e^{−iσ Re ℓ_t}),
  λ(ξ₀) = inf_σ (1 + X/B).   An den Enden = Lemma O10 (λ_L = 0.98988).
Aufruf: python o26_streifen.py modell
"""
import json
import os
import sys

import numpy as np

import o_kern as K

sys.stdout.reconfigure(encoding="utf-8")
name = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] in ("kneser", "kubisch", "radial", "moeb_ell") else "kneser"
L = K.Lc
Lb = L.conjugate()
R0 = abs(L)
thL = np.angle(L)


def fkt(name):
    if name == "kneser":
        d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "gemeinsam",
                                        "kneser_halbexp_koeff_4000.json"), encoding="utf-8"))
        a = np.array([float(x) for x in d["a"]])
        k = np.arange(len(a))
        f = lambda x: np.polynomial.polynomial.polyval(x, a)
        fp = lambda x: np.polynomial.polynomial.polyval(x, a[1:] * k[1:])
        return f, fp, (lambda x: np.exp(x) / fp(x))
    f, fp = K.f_numpy(name)
    return f, fp, (lambda x: fp(f(x)))


f, fp, mfun = fkt(name)
M = lambda x: (x - L) / (x - Lb)
Mp = lambda x: (L - Lb) / (x - Lb) ** 2
th0 = np.angle(M(-R0))          # Γin enthält −ρ


def zeta(x):
    z = np.log(M(x)) - 1j * th0
    return z


def Minv(w):                       # w = e^{ζ + iθ₀}
    return (L - Lb * w) / (1 - w)


def daten(xi, oben):
    eta = np.pi if oben else 0.0
    zz = xi + 1j * eta
    x = Minv(np.exp(zz + 1j * th0))
    fx = f(x)
    gz = np.log(M(fx)) - 1j * th0
    # Zweig: g(ζ) − ζ nahe ℓ (Im zwischen 0 und π)
    d = gz - zz
    d = d.real + 1j * ((d.imag + np.pi) % (2 * np.pi) - np.pi)
    gp = (Mp(fx) * fp(x) / M(fx)) * (M(x) / Mp(x))
    m = mfun(x)
    Rw = np.sqrt(gp / fp(x)) / m
    # Zweig von sqrt: stetig, Start am L-Ende bei μ^{-3/2}
    return x, d, Rw, np.abs(m) ** 2


mu = np.sqrt(L)
ell = np.log(mu)


def haupt():
    xis = np.linspace(-14, 14, 561)
    mu = np.sqrt(L)
    ell = np.log(mu)
    print(f"{name}: ℓ = log μ = {ell:.6f},  μ^(-3/2) = {mu ** -1.5:.6f}")
    rows = []
    for xi in xis:
        x_b, lb, Rb, pb = daten(xi, False)
        x_t, lt, Rt, pt = daten(xi, True)
        rows.append((xi, lb, lt, Rb, Rt, pb, pt))
    # Zweigwahl für R stetig machen (Vorzeichen der Wurzel)
    for j in (3, 4):
        prev = mu ** -1.5
        for i in range(len(rows)):
            r = list(rows[i])
            if abs(r[j] - prev) > abs(-r[j] - prev):
                r[j] = -r[j]
            prev = r[j]
            rows[i] = tuple(r)
    sig = np.linspace(-40, 40, 16001)


    def lam_frozen(lb, lt, Rb, Rt, pb, pt):
        B = pb + pt * abs(Rt) ** 2 * np.exp(-2 * sig * (np.pi + lt.imag))
        X = pb * np.real(Rb * np.exp(1j * sig * lb)) + pt * np.exp(-sig * (2 * np.pi + lt.imag)) * np.real(np.conj(Rt) * np.exp(-1j * sig * lt.real))
        q = 1 + X / B
        i = np.argmin(q)
        return q[i], sig[i], q.max()


    print(" ξ₀      ℓ_b                ℓ_t                R_b                 R_t                p_b     p_t     λ(ξ₀)    σ*      sup q")
    lams = []
    for (xi, lb, lt, Rb, Rt, pb, pt) in rows:
        lmin, smin, qmax = lam_frozen(lb, lt, Rb, Rt, pb, pt)
        lams.append(lmin)
        if abs(xi * 2 - round(xi * 2)) < 1e-9:
            print(f"{xi:+6.2f}  {lb.real:+.4f}{lb.imag:+.4f}i  {lt.real:+.4f}{lt.imag:+.4f}i  {Rb.real:+.4f}{Rb.imag:+.4f}i  "
                  f"{Rt.real:+.4f}{Rt.imag:+.4f}i  {pb:.4f}  {pt:.4f}  {lmin:.5f}  {smin:+6.2f}  {qmax:.3f}")
    print(f"min_ξ₀ λ(ξ₀) = {min(lams):.6f} bei ξ₀ = {xis[int(np.argmin(lams))]:+.2f}")
    # Ableitungen (Variation der Koeffizienten)
    arr = np.array([[r[1], r[2], r[3], r[4]] for r in rows])
    dx = xis[1] - xis[0]
    for j, nm in enumerate(["ℓ_b", "ℓ_t", "R_b", "R_t"]):
        der = np.abs(np.diff(arr[:, j])) / dx
        print(f"max |d{nm}/dξ| = {der.max():.4f} bei ξ = {xis[int(np.argmax(der))]:+.2f}")


if __name__ == "__main__":
    haupt()

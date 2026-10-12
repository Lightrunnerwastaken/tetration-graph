"""k53 — Test der Dzjadyk-Rate an einer einspringenden Ecke (BERICHT §35), unabhängig von f_K:
g(z) = (z − L)^{s0} (Zweig mit Schnitt in der äusseren Lücke), E = Kompaktum mit Ecke in L, Aussenöffnung Θ (zwei Strahlen θ_a, θ_b), sonst glatter Rand.
Polynom-Bestapproximation näherungsweise per kleinsten Quadraten auf dichten Randpunkten (Vandermonde mit Arnoldi, Brubeck–Nakatsukasa–Trefethen),
gemessen sup-Fehler auf feineren Randpunkten; Vorhersage dist ~ N^{−Re(s0)·Θ/π}.
Aufruf: python k53_dzjadyk_test.py
"""
import math
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")


def rand(L, th_a, th_b, R, m):
    """Randpunkte: Strahl θ_b von L nach aussen (Länge R), Bogen, Strahl θ_a zurück nach L; Innen = Seite mit Öffnung 360° − Θ."""
    # Öffnung des Innern: von θ_a gegen den Uhrzeigersinn bis θ_b + 360 (θ_b < θ_a)
    t = np.linspace(0, 1, m)
    ray_b = L + np.outer(t ** 2, [1])[:, 0] * R * np.exp(1j * th_b)      # verdichtet zur Ecke
    phis = np.linspace(th_b, th_a - 2 * math.pi, 3 * m)                  # Bogen im Uhrzeigersinn durch das Innere-Gegenüber
    # Bogen: Kreis um L mit Radius R, Winkel von θ_b abwärts bis θ_a − 2π (also durch die Innenöffnung)
    arc = L + R * np.exp(1j * phis)
    ray_a = L + (t[::-1] ** 2) * R * np.exp(1j * th_a)
    return np.concatenate([ray_b, arc, ray_a])


def g(z, L, s0, th_cut):
    w = z - L
    ang = np.angle(w * np.exp(-1j * th_cut))           # Schnitt bei θ_cut: Winkel in (−π, π] relativ dazu
    ang = np.where(ang <= 0, ang + 2 * math.pi, ang)    # Zweig: (0, 2π] ab dem Schnitt
    r = np.abs(w)
    lw = np.log(np.maximum(r, 1e-300)) + 1j * (ang + th_cut)
    return np.where(r > 0, np.exp(s0 * lw), 0)


def va_fit(z, f, n, zc, scale):
    """Vandermonde mit Arnoldi: Basis q_k orthonormal auf den Stützstellen."""
    x = (z - zc) / scale
    M = len(z)
    Q = np.zeros((M, n + 1), dtype=complex); H = np.zeros((n + 1, n), dtype=complex)
    Q[:, 0] = 1
    for k in range(n):
        v = x * Q[:, k]
        for j in range(k + 1):
            H[j, k] = np.vdot(Q[:, j], v) / M
            v = v - H[j, k] * Q[:, j]
        H[k + 1, k] = np.linalg.norm(v) / math.sqrt(M)
        Q[:, k + 1] = v / H[k + 1, k]
    c = np.linalg.lstsq(Q, f, rcond=None)[0]
    return c, H


def va_eval(c, H, z, zc, scale):
    x = (z - zc) / scale
    n = H.shape[1]
    W = np.zeros((len(z), n + 1), dtype=complex); W[:, 0] = 1
    for k in range(n):
        v = x * W[:, k]
        for j in range(k + 1):
            v = v - H[j, k] * W[:, j]
        W[:, k + 1] = v / H[k + 1, k]
    return W @ c


def main():
    L = 0.31813150520476413 + 1.3372357014306895j
    s0 = 5.446950720067007 + 1.057939991156939j
    R = 1.6
    for Theta_deg in (111.7, 131.7, 166.7):
        th_b = math.radians(20.0)
        th_a = th_b + math.radians(Theta_deg)           # Aussenlücke (θ_b, θ_a)
        th_cut = 0.5 * (th_a + th_b)
        zc = L + 0.6 * R * np.exp(1j * (th_a + math.pi - math.radians(Theta_deg) / 2 + math.pi))  # grob im Innern
        zc = L - 0.5 * R * np.exp(1j * th_cut)
        zs = rand(L, th_a, th_b, R, 1500)
        zt = rand(L, th_a, th_b, R, 6000)
        fs = g(zs, L, s0, th_cut); ft = g(zt, L, s0, th_cut)
        rows = []
        prev = None
        for n in (10, 20, 30, 40, 60, 80):
            c, H = va_fit(zs, fs, n, zc, R)
            err = np.max(np.abs(va_eval(c, H, zt, zc, R) - ft))
            loc = "" if prev is None else " (lok. Exp. %.2f)" % (math.log(prev[1] / err) / math.log(n / prev[0]))
            rows.append("N=%d: %.2e%s" % (n, err, loc))
            prev = (n, err)
        alpha = Theta_deg / 180.0
        print("Θ = %.1f° (α = %.3f), Vorhersage Exponent %.2f: %s" % (Theta_deg, alpha, 5.417 * alpha, " | ".join(rows)), flush=True)


if __name__ == "__main__":
    main()

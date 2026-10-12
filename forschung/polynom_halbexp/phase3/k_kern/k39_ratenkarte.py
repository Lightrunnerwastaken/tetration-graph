"""k39 — Äussere Rate des adjungierten Kopfes bei L aus dem expliziten Eckprofil (Sattel k = −1, der den Wachstumssektor an der Γin-Seite trägt, k34)
und Lage der Null-Spirale relativ zur Kreistangente und zu Keil-Linsen.

Mit u = Λ₀ n^{a}(x − L), Λ₀ = −λ_Aρ/L² (Gs Koordinate, linearisiert), ist log|F_A(u)| ≈ Re(−a s*) + Re((β_A − ½) log s*), s* = (u (iα)^{−iα} e^{−2πi})^{κ₁}.
Re(−a s*)/n = g(x) ist n-unabhängig (Kegelfunktion Re(C(x−L)^{κ₁})): die n-te-Wurzel-Rate des Kopfes relativ zu seinem Eckwert.
Ausgabe: (1) Kontrolle gegen die Reihe (|u| = 40, 80); (2) g längs Richtungen θ = θ_c − φ (θ_c = arg L + 90° = Γin-Tangente, φ = Keilwinkel nach aussen)
für r = 10^{−1} … 10^{−40}; Nullstellen r_×(φ); N_×(φ) := (ρ/r_×)^{1/Re a} (Kegelskala erreicht r_×); (3) Drehrate der Null-Spirale.
Aufruf: python k39_ratenkarte.py
"""
import sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 60
L = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237", "1.33723570143068940890116214319371061253950213846051241687")
RHO = abs(L)
ALPHA = L / (4 * mp.pi)
A = 1 + 1j * ALPHA
KAP = 1 / A
BA = 0.5 - 2j * ALPHA
LAM = mp.mpc(-0.203297, -1.406008)
LAM0 = -LAM * RHO / L ** 2
THC = mp.arg(L) + mp.pi / 2


def sstar(u, k=-1):
    return mp.exp(KAP * (mp.log(u) + 2j * mp.pi * k - 1j * ALPHA * mp.log(1j * ALPHA)))


def rate(r, theta, n=mp.mpf(10) ** 8, k=-1):
    x_L = r * mp.expjpi(theta / mp.pi)
    u = LAM0 * mp.power(n, A) * x_L
    return mp.re(-A * sstar(u, k)) / n, u


def main():
    K = 900
    cY = [mp.gamma(BA - 1j * ALPHA * kk) / (mp.factorial(kk) * mp.gamma(BA)) for kk in range(K)]

    def F(u):
        s = mp.mpc(0); p = mp.mpc(1)
        for c in cY:
            s += c * p; p *= -u
        return s
    print("(1) Kontrolle Sattel k=−1 gegen Reihe: log|F| − Re(−as*) − Re((β−½)log s*):")
    for ang in [100, 110, 120, 140, 160, 180]:
        row = []
        for R in [40, 80]:
            u = R * mp.expjpi(ang / 180)
            ss = sstar(u)
            row.append("%.3f" % float(mp.log(abs(F(u))) - mp.re(-A * ss) - mp.re((BA - 0.5) * mp.log(ss))))
        print("   arg u = %d°: %s" % (ang, " / ".join(row)))
    print("\n(2) Rate g(r) = Re(−a s*_{−1})/n längs θ = θ_c − φ (θ_c = %.3f°); n-Unabhängigkeit: n = 1e8 und 1e20" % float(mp.degrees(THC)))
    for phi in [0, 5, 10, 15, 20, 26, 30, 38.3, 45, 60, 90]:
        th = THC - mp.radians(phi)
        row = []
        sign_prev = None
        rx = None
        for e in range(1, 61):
            r = mp.mpf(10) ** (-e / 1.5)
            g, u = rate(r, th)
            if sign_prev is not None and (g > 0) != sign_prev and rx is None:
                rx = r
            sign_prev = g > 0
        g8 = [rate(mp.mpf(10) ** (-ee), th)[0] for ee in (1, 2, 3, 5, 10, 20)]
        g20 = rate(mp.mpf(10) ** -3, th, n=mp.mpf(10) ** 20)[0]
        argu = mp.degrees(mp.arg(rate(mp.mpf(10) ** -3, th)[1]))
        nx = (RHO / rx) ** (1 / mp.re(A)) if rx is not None else None
        print("   φ = %5.1f°: g(10^-1,-2,-3,-5,-10,-20) = %s | g(1e-3; n=1e20) = %.3e | arg u(1e-3, n=1e8) = %.1f° | Vorzeichenwechsel bei r_× ≈ %s, N_× ≈ %s"
              % (phi, " ".join("%.2e" % float(v) for v in g8), float(g20), float(argu), mp.nstr(rx, 3) if rx else "—", mp.nstr(nx, 3) if nx else "—"))
    # Drehrate der Nullspirale: Re(−a s*) = 0 ⇔ arg(−a s*) = ±π/2; arg(s*) = Im(κ₁)ln|x−L| + Re(κ₁)θ + const ⇒ dθ/dln r = −Im κ₁/Re κ₁
    print("\n(3) Drehrate der Null-Spiralen: dθ/d ln r = −Im κ₁/Re κ₁ = %.5f rad (= %.3f° je e-Faktor in r; θ fällt, wenn r fällt)"
          % (float(-mp.im(KAP) / mp.re(KAP)), float(mp.degrees(-mp.im(KAP) / mp.re(KAP)))))
    # Position der Nullspirale relativ zur Tangente als Funktion von r
    print("    Null-Spirale (Wurzel von g längs Kreisen |x−L| = r, Suchbereich θ ∈ [θ_c − 150°, θ_c + 30°]): Winkel relativ zur Tangente (positiv = nach aussen):")
    for e in [1, 1.5, 2, 3, 5, 10, 20, 40]:
        r = mp.mpf(10) ** (-e)
        f = lambda t: rate(r, THC - mp.radians(t))[0]
        roots = []
        prev = None
        for t in range(-30, 151, 2):
            v = f(t)
            if prev is not None and (v > 0) != (prev[1] > 0):
                roots.append(mp.findroot(f, (prev[0], t), solver="bisect"))
            prev = (t, v)
        print("      r = 1e-%g: Nullstellen bei φ = %s" % (e, ", ".join("%.2f°" % float(z) for z in roots)))


if __name__ == "__main__":
    main()

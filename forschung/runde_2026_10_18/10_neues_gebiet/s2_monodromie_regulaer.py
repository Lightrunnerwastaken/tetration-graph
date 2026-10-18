"""s2: Der reguläre Erzeuger J_L = L σ/σ' (Basis e) entlang von Wegen fortsetzen.

Satz T5: σ und J_L lassen sich längs jedes Weges in C \\ P fortsetzen, P = {0, 1, e, e^e, ...}.
Fortsetzung hier konstruktiv: Zweige der iterierten Logarithmen werden stetig mitgeführt
(ℓ_k(t) = log ℓ_{k-1}(t) + 2πi m, m so, dass ℓ_k stetig bleibt), danach σ = L^n σ_loc(ℓ_n).
Test der Monodromie um 0: Ist 0 ein echter Verzweigungspunkt, so ist J_L entlang L -> L+2πi
an der Stelle L+2πi NICHT null (bei trivialer Monodromie wäre J_L(L+2πi) = e^{-L} J_L(L) = 0).
Gegenprobe: zwei Schrittweiten / zwei Tiefen n.
"""
import sys
from pathlib import Path
import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib10 as lib


def fortsetzen(B, weg, n, schritte):
    """weg: t -> z(t), t in [0,1]; liefert Liste (t, z, J_L(z)) längs des Weges."""
    L = B.L
    ell = None
    out = []
    for j in range(schritte + 1):
        t = mp.mpf(j) / schritte
        z = weg(t)
        neu = []
        w = z
        for k in range(n):
            lw = mp.log(w)
            if ell is not None:
                # Zweig wählen, der dem vorigen Wert am nächsten liegt
                m = mp.nint(mp.im(ell[k] - lw) / (2 * mp.pi))
                lw = lw + 2j * mp.pi * m
            neu.append(lw)
            w = lw
        ell = neu
        u = ell[-1] - L
        su = mp.polyval(B.S[::-1], u)
        dsu = mp.polyval(B.dS[::-1], u)
        # σ(z) = L^n σ_loc(ℓ_n), σ'(z) = L^n σ_loc'(ℓ_n) ∏ 1/ℓ_{k-1}  (ℓ_0 = z)
        prod = mp.mpc(1)
        prev = z
        for k in range(n):
            prod /= prev
            prev = ell[k]
        J = L * su / (dsu * prod)
        out.append((t, z, J, abs(u)))
    return out


if __name__ == "__main__":
    # Aufruf: python s2_monodromie_regulaer.py [1|2]  (Einstellung 1: dps 30, n 60, 4000 Schritte; 2: dps 36, n 75, 2500 Schritte)
    wahl = sys.argv[1] if len(sys.argv) > 1 else "1"
    einst = {"1": (30, 40, 60, 4000), "2": (36, 46, 75, 2500)}[wahl]
    for dps, K, n, schritte in [einst]:
        mp.mp.dps = dps
        B = lib.Basis(K=K)
        L = B.L
        # Weg 1: senkrecht L -> L + 2πi (die Bildkurve e^z umläuft 0 einmal)
        weg1 = lambda t: L + 2j * mp.pi * t
        r = fortsetzen(B, weg1, n, schritte)
        print(f"dps={dps} n={n} Schritte={schritte}: max|ℓ_n - L| längs Weg = {mp.nstr(max(x[3] for x in r), 3)}")
        print("   J_L(L) =", mp.nstr(r[0][2], 15), "  J_L(L+2πi) [fortgesetzt] =", mp.nstr(r[-1][2], 20))
        print("   Vergleich: e^{-L} J_L(L) = 0;  J_L(L+2πi) / L =", mp.nstr(r[-1][2] / L, 20))
        # Weg 2: Kreis um L+2πi... als Kontrolle: geschlossene Schleife um 2πi*0.5 (kein Punkt von P umschlossen
        # im Bild?) -> hier einfache Kontrolle: Schleife um L selbst mit Radius 0.3 (J_L dort holomorph)
        weg2 = lambda t: L + mp.mpf("0.3") * mp.expj(2 * mp.pi * t)
        r2 = fortsetzen(B, weg2, n, schritte)
        print("   Schleife um L (Kontrolle, J_L holomorph): Sprung =", mp.nstr(abs(r2[-1][2] - r2[0][2]), 3))
        # Weg 3: Schleife um 0 mit Radius 1/2, Start bei 1/2 (umschliesst nur 0 aus P)
        weg3 = lambda t: mp.mpf("0.5") * mp.expj(2 * mp.pi * t)
        r3 = fortsetzen(B, weg3, n, schritte)
        print("   Schleife um 0 (r = 1/2): J vorher", mp.nstr(r3[0][2], 15), " nachher", mp.nstr(r3[-1][2], 15))
        # Weg 4: Schleife um 1 mit Radius 1/2 (umschliesst nur 1 aus P)
        weg4 = lambda t: 1 + mp.mpf("0.5") * mp.expj(2 * mp.pi * t)
        r4 = fortsetzen(B, weg4, n, schritte)
        print("   Schleife um 1 (r = 1/2): J vorher", mp.nstr(r4[0][2], 15), " nachher", mp.nstr(r4[-1][2], 15))
        # Weg 5: Schleife um 3 mit Radius 1/4 (umschliesst keinen Punkt aus P) -> Sprung soll 0 sein
        weg5 = lambda t: 3 + mp.mpf("0.25") * mp.expj(2 * mp.pi * t)
        r5 = fortsetzen(B, weg5, n, schritte)
        print("   Schleife um 3 (r = 1/4, kein P-Punkt): Sprung =", mp.nstr(abs(r5[-1][2] - r5[0][2]), 3))
        sys.stdout.flush()

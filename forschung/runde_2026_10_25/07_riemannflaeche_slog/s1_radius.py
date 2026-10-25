"""s1: Taylor-Koeffizienten A_n von Knesers slog um 0 (nur lesend aus polynom_halbexp/gemeinsam).
Prueft: A_n = -2 Re( 1/(n L^(n+1)) ) * (1 + O(n^-Re(gamma)))  -> Radius genau |L|, Log-Koeffizient 1/L bei L.
Ausserdem: slog(w) = A(e^w) - 1 auf Re w < Re L (2 pi i-periodisch), Stichprobe gegen die Taylorreihe selbst."""
import json
import mpmath as mp

mp.mp.dps = 50
D = json.load(open(r"D:\Users\janis\Documents\Matura\tetration-graph\forschung\polynom_halbexp\gemeinsam\kneser_slog_koeff_300.json"))
A = [mp.mpf(a) for a in D["A"]]
L = mp.mpc("0.3181315052047641353126542515876645172035", "1.33723570143068940890116214319371061254")
g = 2j * mp.pi / L
print("|L| =", mp.nstr(abs(L), 20), " e^ReL =", mp.nstr(mp.e ** L.real, 20))
print(" n   A_n / (-2Re(1/(n L^(n+1))))   (Abweichung vom Leitglied relativ zu |L|^-n/n)")
for n in [20, 40, 60, 80, 100, 150, 200, 250, 300]:
    lead = -2 * mp.re(1 / (n * L ** (n + 1)))
    rel = (A[n] - lead) / (abs(L) ** (-n - 1) / n)
    print(f"{n:4d}  ratio={mp.nstr(A[n]/lead, 12):>16}  rest*n|L|^(n+1)={mp.nstr(abs(rel), 6):>10}  n^-Re(g)={mp.nstr(mp.mpf(n)**(-g.real),6)}")
# Fit: |rest| ~ C n^-p
import math
pts = []
for n in range(100, 301, 20):
    lead = -2 * mp.re(1 / (n * L ** (n + 1)))
    rel = abs((A[n] - lead) / (abs(L) ** (-n - 1) / n))
    pts.append((math.log(n), float(mp.log(rel))))
# Huellkurve: Maximum je Fenster (Rest oszilliert)
print("Rest-Exponent (lokale Steigung der Huellkurve):")
env = []
for n0 in range(60, 300, 40):
    m = max(abs((A[n] - (-2 * mp.re(1 / (n * L ** (n + 1))))) / (abs(L) ** (-n - 1) / n)) for n in range(n0, n0 + 40))
    env.append((n0 + 20, m))
for (n1, m1), (n2, m2) in zip(env, env[1:]):
    print(f"  n~{n1}-{n2}: Steigung {mp.nstr(mp.log(m2/m1)/mp.log(mp.mpf(n2)/n1), 6)}  (Theorie -Re(gamma) = {mp.nstr(-g.real, 6)})")


def Aser(z, N=300):
    return mp.polyval(A[N::-1], z)


# Periodische Halbebene: slog(w) = A(e^w) - 1, Stichprobe: Konsistenz mit Taylorreihe selbst in |w|<0.3
for w in [mp.mpc(-0.2, 0.15), mp.mpc(0.1, -0.25)]:
    print("Konsistenz A(w) vs A(e^w)-1 :", mp.nstr(abs(Aser(w) - (Aser(mp.exp(w)) - 1)), 5))

"""Runde 2, Schritt 4: Welches Gesetz regiert die Pol-Antwort im Möbius-Modell — Cs Max-Formel max(1/|ζ|, R₀/|ζ−a₀|) oder
die Green-Funktion von ℂ∖Ḡ' (Bernstein–Walsh)?  Green-Funktion explizit über Cayley + Potenz: Ḡ' ↦ Sektor, Komplement Σ ↦ H."""
import sys, time
import mpmath as mp
from flint import acb
from moebius_lin import setup, pol_antwort, auswerten

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 30
c = mp.mpf(sys.argv[1])
Ns = [int(x) for x in sys.argv[2].split(",")]
tau = mp.tan(c / 2)
R0 = mp.sqrt(1 + tau ** 2)
om = mp.pi - c / 2


def green(zeta):
    W = (1 + 1j * zeta) / (1 - 1j * zeta)
    a = mp.arg(W * mp.expj(-(mp.pi / 2 + c / 2)))
    if a < 0:
        a += 2 * mp.pi
    if not (0 < a < om):
        return mp.mpf(0)
    u = abs(W) ** (mp.pi / om) * mp.expj(a * mp.pi / om)
    a_inf = (mp.pi / 2 - c / 2)          # Argument von −1 nach Drehung
    u_inf = mp.expj(a_inf * mp.pi / om)
    return mp.log(abs((u - mp.conj(u_inf)) / (u - u_inf)))


K = 256
zs = [acb(mp.nstr(mp.cos(2 * mp.pi * j / K), 25), mp.nstr(mp.sin(2 * mp.pi * j / K), 25)) for j in range(K)]
punkte = [mp.mpc(1 / tau, 0), mp.mpc(-2, 0), mp.mpc(2.0, 0), mp.mpc(0, 1.6), mp.mpc(1.2, 1.2), mp.mpc(-1.3, 0.8)]
print(f"c = {c}: τ = {mp.nstr(tau, 5)}, R₀ = {mp.nstr(R0, 5)}")
for zeta in punkte:
    g = green(zeta)
    mx = max(1 / abs(zeta), R0 / abs(zeta - tau))
    print(f"ζ = {mp.nstr(zeta, 4)}: e^(−g) = {mp.nstr(mp.exp(-g), 6)}, Max-Formel = {mp.nstr(mx, 6)}")
res = {}
for N in Ns:
    t0 = time.time()
    S = setup(c, N, 900)
    Ji = S["J"].inv()
    for zeta in punkte:
        za = acb(mp.nstr(mp.re(zeta), 25), mp.nstr(mp.im(zeta), 25))
        coef = pol_antwort(S, za, Ji)
        best = mp.mpf(0)
        for z in zs:
            v = auswerten(coef, z) - 1 / (za - z)
            best = max(best, abs(mp.mpc(v.real.mid().str(20, radius=False), v.imag.mid().str(20, radius=False))))
        res[(N, zeta)] = best
    print(f"N = {N} fertig ({time.time() - t0:.0f}s)", flush=True)
for zeta in punkte:
    g = green(zeta)
    mx = max(1 / abs(zeta), R0 / abs(zeta - tau))
    zeile = []
    for N in Ns:
        zeile.append(f"N={N}: |E|={mp.nstr(res[(N, zeta)], 3)} /e^(−Ng)={mp.nstr(res[(N, zeta)] * mp.exp(N * g), 3)} /max^N={mp.nstr(res[(N, zeta)] / mx ** N, 3)}")
    print(f"ζ = {mp.nstr(zeta, 4)}: " + "  ".join(zeile))

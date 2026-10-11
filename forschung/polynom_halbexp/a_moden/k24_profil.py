"""Strang A, Runde 2: Eckschicht-Profil von E_N = Π_N f_K − f_K an L.

E_N(L(1 − v/N^α)) für α = 1 und α = 0.88, v = ρ e^{iφ} (ρ-Gitter, φ ∈ [−90°, 90°]; φ = ±90° ist der Rand |z| ≈ |L|),
für N aus daten/lin_PN.json. Normierung Φ_N(v) = E_N · N^{αγ} (Ansatz E ≈ ℓ^γ Φ(v), ℓ = N^{−α}).
Ausgabe: gemeinsam/eckschicht_profil.json (Rohwerte E_N und Φ_N), Konvergenztabelle im Log.
"""
import json
import sys
import time

sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402
from flint import acb, acb_poly, arb, ctx  # noqa: E402

from lin_kern import DATEN, GEM, lade  # noqa: E402

mp.mp.dps = 45
ctx.prec = 200
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
G = 2j * mp.pi / L
SL = mp.exp(L / 2)
LML = mp.log(-L)
cx = lambda p: mp.mpc(mp.mpf(p[0]), mp.mpf(p[1]))  # noqa: E731
sth = json.loads((DATEN / "singularitaet_theta.json").read_text(encoding="utf-8"))
phi = [cx(p) for p in sth["phi"]]
psi = [cx(p) for p in sth["psi"]]
eps = [cx(p) for p in sth["eps_q"]]


def pv(c, x):
    s = mp.mpc(0)
    for a in reversed(c):
        s = s * x + a
    return s


def fK_lokal(z):
    t = z - L
    if t == 0:
        return L
    s = pv(psi, t)
    logs = LML + mp.log(1 - z / L) + mp.log(s / t)
    qq = mp.exp(G * logs)
    e = mp.fsum(eps[k] * qq ** k for k in range(1, len(eps)))
    return pv(phi, SL * s * mp.exp(L * e))


RHOS = [0, 0.1, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32]
PHIS = [-90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90]
ALPHAS = [1.0, 0.88]
Ldat = lade()
Ns = [N for N in sorted(Ldat) if N in (100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000)]
pfad = GEM / "eckschicht_profil.json"
aus = json.loads(pfad.read_text(encoding="utf-8")) if pfad.is_file() else {
    "beschreibung": ("Eckschicht von E_N = Π_N f_K − f_K an L (Strang A, Runde 2). Punkte z = L(1 − v/N^α), v = ρ e^{iφ} "
                     "(φ in Grad; φ = ±90° ≈ Rand |z| = |L|, φ = 0 radial nach innen). Pro α, N, ρ, φ: E = [Re, Im] und "
                     "Φ = E·N^{αγ}, γ = 2πi/L. f_K nahe L aus der exakten lokalen Formel (k5/k6), P_N = Π_N f_K aus k21."),
    "gamma": [mp.nstr(mp.re(G), 20), mp.nstr(mp.im(G), 20)], "L": [mp.nstr(mp.re(L), 20), mp.nstr(mp.im(L), 20)],
    "daten": {}}
t0 = time.time()
for N in Ns:
    for al in ALPHAS:
        key = f"{al}|{N}"
        if key in aus["daten"]:
            continue
        bp = acb_poly([acb(arb(x)) for x in Ldat[N]["b"]])
        ell = mp.power(N, -al)
        norm = mp.power(N, al * G)
        tab = {}
        for rho in RHOS:
            for ph in (PHIS if rho > 0 else [0]):
                v = mp.mpf(rho) * mp.expj(mp.pi * ph / 180)
                z = L * (1 - v * ell)
                if abs(z - L) > mp.mpf("0.12"):
                    continue
                w = bp(acb(arb(mp.nstr(mp.re(z), 42)), arb(mp.nstr(mp.im(z), 42))))
                E = mp.mpc(mp.mpf(w.real.mid().str(40, radius=False)), mp.mpf(w.imag.mid().str(40, radius=False))) - fK_lokal(z)
                Ph = E * norm
                tab[f"{rho}|{ph}"] = {"E": [mp.nstr(mp.re(E), 12), mp.nstr(mp.im(E), 12)], "Phi": [mp.nstr(mp.re(Ph), 10), mp.nstr(mp.im(Ph), 10)]}
        aus["daten"][key] = tab
        pfad.write_text(json.dumps(aus, indent=0, ensure_ascii=False), encoding="utf-8")
        print(f"α = {al}, N = {N}: {len(tab)} Punkte ({time.time() - t0:.0f} s)", flush=True)
# Konvergenztabelle radial (φ = 0) und am Rand (φ = −90°)
for al in ALPHAS:
    print(f"\nα = {al}: Φ_N(ρ, φ) = E·N^(αγ), Spalten N")
    for rho, ph in [(0, 0), (0.5, 0), (1, 0), (2, 0), (4, 0), (8, 0), (16, 0), (1, -90), (2, -90), (1, 90), (4, 45), (4, -45)]:
        zeile = []
        for N in Ns:
            d = aus["daten"].get(f"{al}|{N}", {}).get(f"{rho}|{ph}")
            zeile.append(mp.nstr(cx(d["Phi"]), 4) if d else "—")
        print(f"  ρ={rho:>4}, φ={ph:>4}: " + "  ".join(f"{s:>24}" for s in zeile))

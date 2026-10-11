"""Strang P, Satz P-V Schritt 3: RIGOROSES Zertifikat für Knesers sexp (Basis e) um die Daten (T̃, θ̃).

Geometrie (Ebene): Schleife Γ₊ = [−½, ½] + iH, Naht Σ = −½ + i[−H, H], Γ₋ = Spiegel; Nachbarschaft ρ.
  G_B = [−0.6−g, 0.6+g] × [H−ρ−g, H+ρ+g]  (Band),  G_S = [−0.5−ρ−g, −0.5+ρ+g] × [−H−ρ−g, H+ρ+g]  (Naht), g = kleiner Rand.
Schritte (alles in arb):
 (1) η_B ≥ sup_{G_B}|Ũ − T̃|, η_S ≥ sup_{G_S}|T̃(·+1) − e^{T̃}|  (Maximumprinzip: Taylor-Modelle auf dem Rand der Rechtecke).
 (2) |Ũ'| ≥ m_U, |Ũ''| ≤ M_U auf G_B ⊕ D(s₀); |T̃'| ≥ m_T, |T̃''| ≤ M_T auf (G_S + 1) ⊕ D(s₀)  (grobe Ball-Überdeckung).
     ⇒ e (Ũ(ζ+e) = T̃(ζ)) und d (T̃(ζ+1+d) = e^{T̃(ζ)}) existieren auf G_B bzw. G_S, |e| ≤ 2η_B/m_U, |d| ≤ 2η_S/m_T (Kantorovich skalar).
 (3) Netz-Konstante K₀ (Zylinder-Cauchy-Integrale über Γ±, Σ mit knotenverträglichen Dichten, |Dichte| ≤ 1 auf der ρ-Nachbarschaft).
 (4) Fixpunkt-Abbildung W ↦ Cousin(e(·+W), d(·+W), Rampe) ist Kontraktion auf {‖W‖ ≤ r}; Radius und Lipschitz-Konstante.
Ausgabe: p27_zert.json (alle Schranken). Die Identifikation (Paulsen–Cowgill) prüft p28.
Aufruf: python p27_kneser_zert.py DATEN.json [prec] [m]
"""
import json
import sys
import time
from pathlib import Path

from flint import acb, arb, ctx

import kc_rig

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
DATEN = HIER / sys.argv[1]
PREC = int(sys.argv[2]) if len(sys.argv) > 2 else 900
M = int(sys.argv[3]) if len(sys.argv) > 3 else 240
H = arb("0.40")
RHO = arb("0.10")
G = arb("0.01")
RT = arb("0.025")          # Taylor-Modell-Radius
RC = arb("0.15")           # Cauchy-Radius für den Schwanz
S0 = arb("0.02")           # Rouché-Radius


def amax(a, b):
    return (a + b + abs(a - b)) / 2


def rand_punkte(x0, x1, y0, y1, schritt):
    """Mittelpunkte auf dem Rechteckrand, Abstand ≤ schritt (Scheiben vom Radius RT überdecken den Rand)."""
    pts = []
    for (za, zb) in [((x0, y0), (x1, y0)), ((x1, y0), (x1, y1)), ((x1, y1), (x0, y1)), ((x0, y1), (x0, y0))]:
        a = acb(arb(za[0]), arb(za[1]))
        b = acb(arb(zb[0]), arb(zb[1]))
        laenge = float(abs(b - a).mid())
        n = max(1, int(laenge / schritt) + 1)
        for j in range(n):
            pts.append(a + (b - a) * j / n)
    return pts


def main():
    t_start = time.time()
    R = kc_rig.Rig(PREC, n0=400, m_extra=300)
    print(f"Rig: C = {R.C.str(4)}, Lipschitz {R.lip.str(3)}  ({time.time() - t_start:.0f}s)")
    d = json.loads(DATEN.read_text(encoding="utf-8"))
    t = [acb(arb(arb(a).mid())) for a, b in d["t"]]
    c = [acb(arb(arb(a).mid()), arb(arb(b).mid())) for a, b in d["c"]]
    raus = {"daten": DATEN.name, "prec": PREC, "m": M, "H": "0.40", "rho": "0.10", "G_B": [-0.55, 0.55, 0.30, 0.50], "G_S": [-0.65, -0.45, -0.55, 0.55]}
    # ---------------- (1) η_B, η_S
    x0, x1 = -0.55, 0.55                       # G_B (Band für e): Re ∈ [−½−0.05, ½+0.05], Im ∈ [H−ρ/2−0.05, H+ρ/2+0.05]
    y0, y1 = float((H - RHO / 2 - arb("0.05")).mid()), float((H + RHO / 2 + arb("0.05")).mid())
    schritt = float((RT * arb("1.8")).mid())
    etaB = arb(0)
    for j, zc in enumerate(rand_punkte(x0, x1, y0, y1, schritt)):
        a, tau = R.U_taylor(c, zc, M, RT, RC)
        b = R.T_taylor(t, zc)
        s = sum((abs(a[i] - b[i]) * RT ** i for i in range(M + 1)), arb(0))
        s += sum((abs(b[i]) * RT ** i for i in range(M + 1, len(b))), arb(0))
        etaB = amax(etaB, ob := (s + tau).abs_upper())
        if j % 10 == 0:
            print(f"  G_B Rand {j}: {complex(zc.real.mid(), zc.imag.mid())}: {float(ob.mid()):.3e}  ({time.time() - t_start:.0f}s)")
            sys.stdout.flush()
    print(f"η_B ≤ {etaB.str(4)}")
    raus["eta_B"] = etaB.str(6)
    xs0, xs1 = float((arb("-0.5") - RHO - arb("0.05")).mid()), -0.45   # G_S (Naht für d)
    ys1 = float((H + RHO / 2 + arb("0.05")).mid())
    etaS = arb(0)
    for j, zc in enumerate(rand_punkte(xs0, xs1, 0.0, ys1, schritt)):     # obere Hälfte; untere per Spiegelung
        a, tau = R.expT_taylor(t, zc, M, RT, RC)
        b = R.T_taylor(t, zc + 1)
        s = sum((abs(b[i] - a[i]) * RT ** i for i in range(M + 1)), arb(0))
        s += sum((abs(b[i]) * RT ** i for i in range(M + 1, len(b))), arb(0))
        etaS = amax(etaS, ob := (s + tau).abs_upper())
        if j % 10 == 0:
            print(f"  G_S Rand {j}: {complex(zc.real.mid(), zc.imag.mid())}: {float(ob.mid()):.3e}  ({time.time() - t_start:.0f}s)")
            sys.stdout.flush()
    print(f"η_S ≤ {etaS.str(4)}  (Rand von G_S: obere Hälfte + Spiegel)")
    raus["eta_S"] = etaS.str(6)
    (HIER / "p27_zert.json").write_text(json.dumps(raus, indent=1), encoding="utf-8")
    # ---------------- (2) Ableitungsschranken (grobe Überdeckung)
    h = arb("0.02")

    def kaestchen(xa, xb, ya, yb):
        """(Mitte, Radius) einer Überdeckung des Rechtecks durch Scheiben (Radius = halbe Breite + halbe Höhe einer Zelle)."""
        out = []
        nx = int(float(((xb - xa) / h).mid())) + 1
        ny = int(float(((yb - ya) / h).mid())) + 1
        for i in range(nx):
            for j in range(ny):
                cx = xa + (xb - xa) * (2 * i + 1) / (2 * nx)
                cy = ya + (yb - ya) * (2 * j + 1) / (2 * ny)
                out.append((acb(arb(cx.mid()), arb(cy.mid())), (xb - xa) / (2 * nx) + (yb - ya) / (2 * ny) + arb(10) ** -30))
        return out

    # Ũ', Ũ'' über Taylor-Modelle niedriger Ordnung an Zellmitten (Ballarithmetik über Kästchen wickelt zu stark)
    mU, MU2 = None, arb(0)
    for zc, rz in kaestchen(arb(x0) - S0, arb(x1) + S0, arb(y0) - S0, arb(y1) + S0):
        a, tau = R.U_taylor(c, zc, 12, 2 * rz, arb("0.15"))
        var = sum((j * abs(a[j]) * rz ** (j - 1) for j in range(2, 13)), arb(0)) + tau / rz
        lo = abs(a[1]) - var
        mU = lo if mU is None else -amax(-mU, -lo)
        MU2 = amax(MU2, sum((j * (j - 1) * abs(a[j]) * rz ** (j - 2) for j in range(2, 13)), arb(0)) + 2 * tau / rz ** 2)
    print(f"|Ũ'| ≥ {mU.str(4)}, |Ũ''| ≤ {MU2.str(4)} auf G_B ⊕ D(s₀)  ({time.time() - t_start:.0f}s)")
    dt = [n * t[n] for n in range(1, len(t))]
    ddt = [n * dt[n] for n in range(1, len(dt))]
    mT, MT2 = None, arb(0)
    for zc, rz in kaestchen(arb(xs0) + 1 - S0, arb(xs1) + 1 + S0, -arb(ys1) - S0, arb(ys1) + S0):
        Z = kc_rig.disc(zc, rz)
        lo = R.T_ball(dt, Z).abs_lower()
        mT = lo if mT is None else -amax(-mT, -lo)
        MT2 = amax(MT2, R.T_ball(ddt, Z).abs_upper())
    print(f"|T̃'| ≥ {mT.str(4)}, |T̃''| ≤ {MT2.str(4)} auf (G_S+1) ⊕ D(s₀)  ({time.time() - t_start:.0f}s)")
    # Kantorovich skalar: |f(0)| ≤ η, |f'| ≥ m, |f''| ≤ M2 auf |s| ≤ s0; hinreichend 2η/m ≤ s0/2 und 2ηM2/m² ≤ 1/2
    assert 2 * etaB / mU < S0 / 2 and 2 * etaB * MU2 / mU ** 2 < arb("0.5")
    assert 2 * etaS / mT < S0 / 2 and 2 * etaS * MT2 / mT ** 2 < arb("0.5")
    e_max = 2 * etaB / mU
    d_max = 2 * etaS / mT
    print(f"|e| ≤ {e_max.str(4)} auf G_B,  |d| ≤ {d_max.str(4)} auf G_S")
    raus.update(mU=mU.str(6), MU2=MU2.str(6), mT=mT.str(6), MT2=MT2.str(6), e_max=e_max.str(6), d_max=d_max.str(6))
    (HIER / "p27_zert.json").write_text(json.dumps(raus, indent=1), encoding="utf-8")
    print(f"Teil (1)–(2) fertig ({time.time() - t_start:.0f}s)")


if __name__ == "__main__":
    main()

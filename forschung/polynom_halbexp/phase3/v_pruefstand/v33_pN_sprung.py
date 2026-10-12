"""Strang V — hochgenaue p_N für N > 1000 (Auftrag Leitung 23:45): Newton-Sprung von pad(p_{N₁}) auf das System T_{N₂}[P∘P − eˣ] = 0.

Start: p_{N₁} (Datei mit Koeffizienten-Strings), auf Grad N₂ mit Nullen aufgefüllt. Newton mit der Tangente J_{N₂}(P)δ = T[δ∘P + (P'∘P)δ],
alle Reihen exakt modulo x^{N₂+1} (arb_poly, explizite Abschneidung), Mitten der Startkoeffizienten als exakte Zahlen.
Präzision 1.5·N₂ + 400 Bit (die Komposition P∘P löscht bis ~0.36·N₂ Dezimalstellen aus). Abbruch, wenn ‖δ‖_{ℓ¹_ρ} < 1e-100 oder nach 8 Schritten.
Ausgabe: Schrittnormen, Residuum, Datei v33_pN_<N₂>.json mit 400 Stellen je Koeffizient.
Aufruf: python v33_pN_sprung.py QUELLE N1 N2      QUELLE = gemeinsam | v33 (eigene Datei v33_pN_<N1>.json)
"""
import json
import sys
import time
from pathlib import Path

from flint import arb, arb_mat, arb_poly, ctx, acb

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
QUELLE, N1, N2 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
TOL = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-100
bits = int(1.5 * N2) + 400
ctx.prec = bits
rho = abs(acb(-1).lambertw())
if QUELLE == "gemeinsam":
    p_str = json.loads((HIER.parents[1] / "gemeinsam" / "pN_N100_bis_1000.json").read_text(encoding="utf-8"))["p"][str(N1)]
    if not Path(HIER.parents[1] / "gemeinsam" / "pN_N100_bis_1000.json").exists():
        raise SystemExit("Quelle fehlt")
else:
    p_str = json.loads((HIER / f"v33_pN_{N1}.json").read_text(encoding="utf-8"))["p"]
p0 = [arb(arb(s).mid()) for s in p_str] + [arb(0)] * (N2 - N1)
n = N2
ex = [arb(1)]
for k in range(1, n + 1):
    ex.append(ex[-1] / k)


def tr(p):
    return p.truncate(n + 1)


def koeff(p, m=n):
    return [p[i] if i <= p.degree() else arb(0) for i in range(m + 1)]


def komp(qc, P):
    r = arb_poly([0])
    for c in reversed(qc):
        r = tr(r * P) + c
    return r


def l1(v):
    return sum(abs(float(v[i].mid())) * float(rho.mid()) ** i for i in range(len(v)))


P = arb_poly(p0)
t00 = time.time()
for it in range(8):
    t0 = time.time()
    PP = komp(koeff(P), P)
    R = [PP[i] - ex[i] if i <= PP.degree() else -ex[i] for i in range(n + 1)]
    rn = l1(R)
    mP = komp(koeff(P.derivative(), n - 1), P)
    J = arb_mat(n + 1, n + 1)
    pw = arb_poly([1])
    mc = koeff(mP)
    for j in range(n + 1):
        if j > 0:
            pw = tr(pw * P)
        c = koeff(pw)
        for i in range(n + 1):
            J[i, j] = c[i] + (mc[i - j] if i >= j else 0)
    t1 = time.time()
    B = arb_mat(n + 1, 1)
    for i in range(n + 1):
        B[i, 0] = -R[i]
    d = J.solve(B, algorithm="approx")
    dx = [d[i, 0] for i in range(n + 1)]
    nd = l1(dx)
    P = P + arb_poly([arb(x.mid()) for x in dx])
    print(f"N = {N1} → {N2}: Schritt {it + 1}: ‖Residuum‖_ℓ¹ρ = {rn:.3e}, ‖δ‖_ℓ¹ρ = {nd:.3e}  [Aufbau {t1 - t0:.0f} s, Lösen {time.time() - t1:.0f} s, {bits} Bit]",
          flush=True)
    if nd < TOL:
        break
PP = komp(koeff(P), P)
R = [PP[i] - ex[i] if i <= PP.degree() else -ex[i] for i in range(n + 1)]
pc = koeff(P)
diff_alt = l1([pc[i] - (p0[i] if i < len(p0) else 0) for i in range(n + 1)])
print(f"Ende: ‖T_N[P∘P − eˣ]‖_ℓ¹ρ = {l1(R):.3e}; ‖p_{N2} − pad(p_{N1})‖_ℓ¹ρ = {diff_alt:.3e}; Zeit {time.time() - t00:.0f} s", flush=True)
(HIER / f"v33_pN_{N2}.json").write_text(json.dumps({"beschreibung": f"p_{N2} per Newton-Sprung von p_{N1} (V, v33), {bits} Bit",
                                                    "p": [x.mid().str(400, radius=False) for x in pc]}), encoding="utf-8")
print(f"gespeichert v33_pN_{N2}.json", flush=True)

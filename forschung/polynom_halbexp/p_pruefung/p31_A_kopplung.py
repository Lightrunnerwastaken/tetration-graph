"""Strang P, Prüfung 31 (Meilenstein 6, Punkt 1): Startbedingung (A) — p_N gegen Π_N f_K, rigoros (arb).

Π_N u := J_N⁻¹T_N J_∞u,  J_∞u = u∘f + (f'∘f)u,  J_N = T_N J_∞|_{P_N}   (ARCHITEKTUR §1).
Schlüssel: f∘f = eˣ ⇒ f'∘f = eˣ/f' (Ableitung) ⇒
    J_∞f = eˣ + (eˣ/f')·f = eˣ(1 + f/f'),   J_N[k][j] = [x^k] f^j + [x^{k−j}] (eˣ/f').
Beides braucht NUR a₀ … a_{N+1} (keinen Schwanz): die Potenzreihe von f'∘f bei 0 hängt sonst von ALLEN a_m ab (f(0) = a₀ ≠ 0).
Eingaben: zertifizierte a_k (gemeinsam/kneser_halbexp_koeff_zert_640.json, Satz P-V), Zweig-Boxen p_N (p30_boxen/box_N.json,
= Kette p18 + p24, auf aicgolling reproduziert, identische Radien).
Ausgaben (ℓ¹_ρ-Normen, ρ = |L|, sind obere Schranken für sup_{D̄}):
  d1 = ‖p_N − Π_N f‖,   d2 = ‖Π_N f − T_640 f‖ + Σ_{N<k≤640}|a_k|ρ^k,   d3 = ‖p_N − T_640 f‖ (dito);
  sup_{D̄}|p_N − f_K| ≤ d3 + S_640, S_640 := Σ_{k>640}|b_k| (Schwanz, aus (K2), p32).
Aufruf: python p31_A_kopplung.py N1 N2 ... (Standard: 140 150 160 170 180 193)
"""
import json
import sys
import time
from pathlib import Path

from flint import acb, arb, arb_mat, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
KOEFF = HIER.parent / "gemeinsam" / "kneser_halbexp_koeff_zert_640.json"
NS = [int(t) for t in sys.argv[1:]] or [140, 150, 160, 170, 180, 193]


def ob(x):
    return x.abs_upper()


def dekodiere(x):
    m, r, e = x
    z = arb(10) ** int(e)
    return arb(int(m)) * z + arb(0, 1) * (arb(int(r)) * z)


def lade_a():
    d = json.loads(KOEFF.read_text(encoding="utf-8"))
    return [arb(m) + arb(0, 1) * arb(r) for m, r, _ in d["a"]]


def reihe_mul(u, v, n):
    return [sum((u[i] * v[k - i] for i in range(max(0, k - len(v) + 1), min(k, len(u) - 1) + 1)), arb(0))
            for k in range(n + 1)]


def reihe_inv(d, n):
    r = [1 / d[0]]
    for k in range(1, n + 1):
        s = sum((d[j] * r[k - j] for j in range(1, min(k, len(d) - 1) + 1)), arb(0))
        r.append(-s / d[0])
    return r


def pi_N(a, N):
    """Π_N f_K als Kugelvektor (Länge N+1)."""
    f = a[:N + 1]
    fs = [(k + 1) * a[k + 1] for k in range(N + 1)]          # f' bis x^N (braucht a_{N+1})
    ex = [arb(1)]
    for k in range(1, N + 1):
        ex.append(ex[-1] / k)
    inv_fs = reihe_inv(fs, N)
    fpf = reihe_mul(ex, inv_fs, N)                           # eˣ/f' = f'∘f
    q = reihe_mul(f, inv_fs, N)                              # f/f'
    rhs = reihe_mul(ex, [arb(1) + q[0]] + q[1:], N)           # eˣ(1 + f/f')
    # Potenzen f^j mod x^{N+1}
    P = [[arb(1)] + [arb(0)] * N]
    for _ in range(N):
        P.append(reihe_mul(P[-1], f, N))
    J = arb_mat(N + 1, N + 1)
    for k in range(N + 1):
        for j in range(N + 1):
            J[k, j] = P[j][k] + (fpf[k - j] if j <= k else 0)
    x = J.solve(arb_mat([[t] for t in rhs]))
    return [x[i, 0] for i in range(N + 1)]


def main():
    ctx.prec = 1200
    rho = abs(acb(-1).lambertw())
    a = lade_a()
    K = len(a) - 1
    aw = [a[k] * rho ** k for k in range(K + 1)]
    ergebnisse = []
    for N in NS:
        t0 = time.time()
        ctx.prec = max(1200, 4 * N + 600)
        pi = pi_N(a, N)
        boxdatei = HIER / "p30_boxen" / f"box_{N}.json"
        p = [dekodiere(x) for x in json.loads(boxdatei.read_text(encoding="utf-8"))["box"]]
        assert len(p) == N + 1
        d1 = sum((ob(p[k] - pi[k]) * rho ** k for k in range(N + 1)), arb(0))
        rest = sum((ob(aw[k]) for k in range(N + 1, K + 1)), arb(0))
        d2 = sum((ob(pi[k] - a[k]) * rho ** k for k in range(N + 1)), arb(0)) + rest
        d3 = sum((ob(p[k] - a[k]) * rho ** k for k in range(N + 1)), arb(0)) + rest
        rad_pi = max(float(pi[k].rad() * rho ** k) for k in range(N + 1))
        # Probe: Π_N reproduziert P_N (B1.2) — hier nur Plausibilität: Π_N f nahe T_N f
        e = dict(N=N, d1_p_minus_PiNf=float(d1.upper()), d2_PiNf_minus_T640f=float(d2.upper()),
                 d3_p_minus_T640f=float(d3.upper()), rest_N_bis_640=float(rest.upper()),
                 radius_PiNf_rho=rad_pi, prec=ctx.prec, t=round(time.time() - t0, 1))
        print(json.dumps(e, ensure_ascii=False), flush=True)
        ergebnisse.append(e)
    alt = {}
    ziel = HIER / "p31_A.json"
    if ziel.exists():
        alt = {e["N"]: e for e in json.loads(ziel.read_text(encoding="utf-8"))["ergebnisse"]}
    for e in ergebnisse:
        alt[e["N"]] = e
    ziel.write_text(json.dumps({"beschreibung": __doc__.splitlines()[0],
                                "ergebnisse": [alt[n] for n in sorted(alt)]}, indent=1, ensure_ascii=False),
                    encoding="utf-8")


if __name__ == "__main__":
    main()

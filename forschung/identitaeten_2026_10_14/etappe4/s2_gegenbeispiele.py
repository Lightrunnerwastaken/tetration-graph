"""Gegenbeispiele: nicht-Kneser-Tetrationen T = sexp(x + p(x)), die eine (oder mehrere) Mittelwert-Identitäten exakt erfüllen.

p = t·(feste Störung) + Σ ε_i·Mode_i, Moden sin 2πkx bzw. cos 2πkx − 1 (p(0) = 0).
Die ε werden mit Newton (Jacobi-Matrix durch zentrale Differenzen, dann Sehnenverfahren) so bestimmt,
dass die gewählten komplexen Bedingungen (Real- und Imaginärteil) verschwinden.

Aufruf: python s2_gegenbeispiele.py FALL DPS K
  FALL = lin        : Jacobi-Matrizen bei Kneser (Linearisierung, Prüfung J_sin = i·J_cos)
  FALL = a|b|c|d    : einzelne Bedingung, Unbekannte Mode 1, feste Störung t·sin 4πx
  FALL = alle       : a,b,c,d und θ̂(−1..−3) = 0 zugleich, Unbekannte Moden 1..7, feste Störung t·sin 16πx
"""
import json
import sys
import time

import mpmath as mp

import e4_lib as E

fall = sys.argv[1]
dps = int(sys.argv[2]) if len(sys.argv) > 2 else 40
k = int(sys.argv[3]) if len(sys.argv) > 3 else 5
E.init(dps)


def baue(fest, moden, eps):
    """fest: Liste (k, a, b); moden: Liste (k, 'sin'|'cos'); eps: Werte."""
    d = {}
    for kk, a, b in fest:
        d.setdefault(kk, [0, 0])
        d[kk][0] += a; d[kk][1] += b
    for (kk, art), e in zip(moden, eps):
        d.setdefault(kk, [0, 0])
        d[kk][0 if art == "sin" else 1] += e
    return ("trig", [(kk, v[0], v[1]) for kk, v in sorted(d.items())])


def residuum(bed, stoer):
    r = E.kenngroessen(stoer, k, moden=6, momente=4)
    # b2: b ohne seinen Anteil erster Ordnung. b = L²V₂ + L⁴(V₄/12 + V₂²/4) + …, also b ≈ L²·a linear (entartet)
    r["b2"] = (r["b"] - E.L ** 2 * r["a"]) / E.L ** 4
    out = []
    for c in bed:
        out += [mp.re(r[c]), mp.im(r[c])]
    return mp.matrix(out), r


def jacobi(bed, fest, moden, eps, h):
    n = len(moden)
    J = mp.matrix(2 * len(bed), n)
    for i in range(n):
        ep = list(eps); em = list(eps)
        ep[i] += h; em[i] -= h
        fp, _ = residuum(bed, baue(fest, moden, ep))
        fm, _ = residuum(bed, baue(fest, moden, em))
        for j in range(2 * len(bed)):
            J[j, i] = (fp[j] - fm[j]) / (2 * h)
    return J


def newton(bed, fest, moden, maxit=40):
    n = len(moden)
    nbed = len(bed)
    eps = [mp.mpf(0)] * n
    h = mp.mpf(10) ** (-(dps // 3))
    t0 = time.time()
    J = jacobi(bed, fest, moden, eps, h)
    print(f"  Jacobi ({time.time()-t0:.0f}s), Konditionszahl {mp.nstr(mp.cond(J), 4)}", flush=True)
    Ji = J ** -1
    for it in range(maxit):
        f, r = residuum(bed, baue(fest, moden, eps))
        nf = mp.norm(f)
        print(f"  it {it}: |F| = {mp.nstr(nf, 3)}", flush=True)
        if nf < mp.mpf(10) ** (-dps + 6):
            break
        if it in (2, 4, 8) or (nbed > 4 and it < 4 and it > 0):   # Jacobi am aktuellen Punkt erneuern
            J = jacobi(bed, fest, moden, eps, h); Ji = J ** -1
        d = Ji * f
        lam = mp.mpf(1)
        while True:   # gedämpfter Schritt
            neu = [e - lam * d[i] for i, e in enumerate(eps)]
            fn, _ = residuum(bed, baue(fest, moden, neu))
            if mp.norm(fn) < nf or lam < mp.mpf(1) / 8:
                break
            lam /= 2
        if lam < 1:
            print(f"    Dämpfung {lam}", flush=True)
        eps = neu
    return eps, r


def bericht(name, stoer, r):
    print(f"== {name}: p = {stoer}")
    for c in ("a", "b", "b2", "c", "d", "neg1", "neg2", "neg3", "V2", "V3", "V4"):
        print(f"   {c:5s} |.| = {mp.nstr(abs(r[c]), 4)}")
    N = mp.fsum(abs(r[f'neg{m}']) ** 2 for m in range(1, 7))
    print(f"   N_6 = Σ_(m≤6)|θ̂_T(−m)|² = {mp.nstr(N, 6)}")


def ser(stoer):
    return [[kk, mp.nstr(a, dps), mp.nstr(b, dps)] for kk, a, b in stoer[1]]


if fall == "lin":
    print(f"Linearisierung bei Kneser (dps {dps}, k {k}): Spalten d/dε für sin 2πmx und cos 2πmx − 1")
    h = mp.mpf(10) ** (-(dps // 3))
    for m in (1, 2, 3):
        for c in ("a", "b", "c", "d", "V3", "neg1", "neg2"):
            vals = {}
            for art in ("sin", "cos"):
                st_p = baue([], [(m, art)], [h]); st_m = baue([], [(m, art)], [-h])
                rp = E.kenngroessen(st_p, k, 3, 3); rm = E.kenngroessen(st_m, k, 3, 3)
                vals[art] = (rp[c] - rm[c]) / (2 * h)
            q = vals["sin"] / vals["cos"] if abs(vals["cos"]) > mp.mpf(10) ** (-dps // 2) else mp.nan
            print(f"  m={m} {c:5s} dC/dε_cos = {mp.nstr(vals['cos'], 12):>40s}   (dC/dε_sin)/(dC/dε_cos) = {mp.nstr(q, 12)}", flush=True)
    sys.exit()

if fall in ("a", "b", "c", "d"):
    t = mp.mpf("0.03")
    fest = [(2, t, 0)]
    moden = [(1, "sin"), (1, "cos")]
    bed = [fall]
elif fall == "alle":
    # a, c, d und θ̂_T(−1..−3) zugleich (b weggelassen: b = L²·a + Glieder höherer Ordnung, siehe Fall ab)
    t = mp.mpf(sys.argv[4]) if len(sys.argv) > 4 else mp.mpf("0.002")
    fest = [(7, t, 0)]
    moden = [(kk, art) for kk in range(1, 7) for art in ("sin", "cos")]
    bed = ["a", "c", "d", "neg1", "neg2", "neg3"]
elif fall == "ab":
    # a und b zugleich; b2 = (b − L²a)/L⁴ ≈ V₄/12 + V₂²/4 trägt den unabhängigen Teil von b
    t = mp.mpf(sys.argv[4]) if len(sys.argv) > 4 else mp.mpf("0.001")
    fest = [(4, t, 0)]
    moden = [(1, "sin"), (1, "cos"), (3, "sin"), (3, "cos")]
    bed = ["a", "b2"]
else:
    raise SystemExit("FALL?")

t0 = time.time()
print(f"Fall {fall}: Bedingungen {bed}, Unbekannte {moden}, fest {fest}, dps {dps}, k {k}")
eps, r = newton(bed, fest, moden)
stoer = baue(fest, moden, eps)
bericht(f"Lösung ({time.time()-t0:.0f}s)", stoer, r)
# min (1 + p')
mn = min(1 + E.p_wert(stoer, mp.mpf(j) / 400)[1] for j in range(400))
print(f"   min(1+p') = {mp.nstr(mn, 6)}")
out = {"fall": fall, "dps": dps, "k": k, "bed": bed, "p": ser(stoer),
       "werte": {c: [mp.nstr(mp.re(v), 20), mp.nstr(mp.im(v), 20)] for c, v in r.items()}}
(E.HIER / "daten" / f"gegenbeispiel_{fall}_{dps}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

"""k60 — K3-2: Ist der Kern an der Ecke lokal integrabel? (CD1) auf der Eckskala zerlegen (BERICHT §42). Kubisch.
(CD1)  (z − t)K_n(z,t) = GEN − NU + XI,  GEN = ρq̃(z)Λ_n(t) (= q_{n+1}(z)λ_n(t)),  NU = t·ν_n(z),  XI = Ξ_n(z,t).
Punkte auf den Spiralen Σ₁ (u_A = 80°, aussen) und Σ₂ (u_A = −90°, innen) in L:
  (a) Eckskala: r = c·|d_n|, c ∈ {0.5, 1, 2, 5}  → Verhältnisse |XI|/|GEN|, |NU|/|GEN| für Paare (z auf Σ₂, t auf Σ₁) und (z, t beide auf Σ₂)
  (b) feste Punkte: r ∈ {0.003, 0.01, 0.03}     → n-Abhängigkeit von XI, NU und der Diagonale GEN(t,t)
Aufruf: python k60_eckkern.py 150 300 600   (schreibt k60.json)
"""
import json
import math
import os
import sys
import time

import mpmath as mp

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, arb_mat, acb_mat, arb_poly, ctx  # noqa: E402
from k_kern import J_tilde, L_RE, L_IM, modell  # noqa: E402
from k58_cd_identitaet import mpc_of, horner  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
HIER = os.path.dirname(os.path.abspath(__file__))


def spirale(L, uA, r, rho):
    th = math.radians(uA + 71.46) + math.radians(1.623) * math.log(r / float(rho))
    return L + r * mp.mpc(math.cos(th), math.sin(th))


def lauf(n):
    t0 = time.time()
    bits = int(1.3 * n) + 400
    ctx.prec = bits
    mp.mp.dps = int(bits / 3.32) + 10
    nd = mp.mp.dps
    L = mp.mpc(L_RE, L_IM); rho = abs(L)
    dn = float(rho) / n
    Jr = J_tilde("kubisch", n, bits=bits); ctx.prec = bits
    fk, mk = modell("kubisch", n)
    rho_a = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
    f = arb_poly(fk[:4])
    pw = arb_poly([1])
    for _ in range(n + 1):
        pw = arb_poly([pw[k] if k <= pw.degree() else arb(0) for k in range(n + 1)]) * f
    B = arb_mat(n + 1, 2)
    B[0, 0] = arb(1)
    for k in range(n + 1):
        ckk = pw[k] if k <= pw.degree() else arb(0)
        B[k, 1] = ckk * rho_a ** k / rho_a ** (n + 1)
    X = Jr.solve(B, algorithm="approx")
    nu = [mp.mpf(X[k, 0].mid().str(nd, radius=False)) for k in range(n + 1)]
    qt = [-mp.mpf(X[k, 1].mid().str(nd, radius=False)) for k in range(n + 1)] + [mp.mpf(1)]
    en = arb_mat(n + 1, 1); en[n, 0] = arb(1)
    Y = Jr.transpose().solve(en, algorithm="approx")
    lam = [mp.mpf(Y[k, 0].mid().str(nd, radius=False)) for k in range(n + 1)]
    mt = [mp.mpf(mk[k].mid().str(nd, radius=False)) * rho ** k for k in range(n + 1)]
    ff = [mp.mpf(fk[k].mid().str(nd, radius=False)) for k in range(4)]
    phit = [(ff[i] - (1 if i == 1 else 0)) * rho ** i for i in range(4)]
    # Punkte
    eck = {("S1", c): spirale(L, 80.0, c * dn, rho) for c in (0.5, 1.0, 2.0, 5.0)}
    eck.update({("S2", c): spirale(L, -90.0, c * dn, rho) for c in (0.5, 1.0, 2.0, 5.0)})
    fest = {("S1", r): spirale(L, 80.0, r, rho) for r in (0.003, 0.01, 0.03)}
    fest.update({("S2", r): spirale(L, -90.0, r, rho) for r in (0.003, 0.01, 0.03)})
    tpunkte = list(eck.items()) + list(fest.items())
    Jc = acb_mat(Jr)
    E = acb_mat(n + 1, len(tpunkte))
    for j, (_, t) in enumerate(tpunkte):
        for k in range(n + 1):
            v = (rho / t) ** k
            E[k, j] = acb(arb(mp.nstr(mp.re(v), nd)), arb(mp.nstr(mp.im(v), nd)))
    R = Jc.solve(E, algorithm="approx")
    rs = [[mpc_of(R[k, j], nd) for k in range(n + 1)] for j in range(len(tpunkte))]
    G = acb_mat(n + 1, len(tpunkte))
    for j, (_, t) in enumerate(tpunkte):
        et = [(rho / t) ** k for k in range(n + 1)]
        comp = [et[k] - mp.fsum(mt[k - i] * rs[j][i] for i in range(k + 1)) for k in range(n + 1)]
        for k in range(n + 1):
            g = mp.fsum(phit[i] * comp[k - i] for i in range(min(3, k) + 1))
            G[k, j] = acb(arb(mp.nstr(mp.re(g), nd)), arb(mp.nstr(mp.im(g), nd)))
    XI = Jc.solve(G, algorithm="approx")
    xis = [[mpc_of(XI[k, j], nd) for k in range(n + 1)] for j in range(len(tpunkte))]
    tindex = {key: j for j, (key, _) in enumerate(tpunkte)}

    def zerlege(zkey, zpkt, tkey):
        j = tindex[tkey]; t = tpunkte[j][1]; z = zpkt
        K = horner(rs[j], z / rho)
        gen = rho * horner(qt, z / rho) * horner(lam, rho / t)
        nuv = t * horner(nu, z / rho)
        xi = horner(xis[j], z / rho)
        lhs = (z - t) * K
        return {"z": str(zkey), "t": str(tkey), "gen": [float(mp.re(gen)), float(mp.im(gen))], "nu": [float(mp.re(nuv)), float(mp.im(nuv))],
                "xi": [float(mp.re(xi)), float(mp.im(xi))], "lhs": [float(mp.re(lhs)), float(mp.im(lhs))],
                "pruef": float(abs(lhs - (gen - nuv + xi)) / max(abs(lhs), mp.mpf(10) ** -30))}

    res = {"n": n, "eck": [], "fest": [], "diag": []}
    for c in (0.5, 1.0, 2.0, 5.0):
        for c2 in (0.5, 1.0, 2.0, 5.0):
            res["eck"].append(zerlege(("S2", c), eck[("S2", c)], ("S1", c2)))
        for c2 in (0.5, 2.0, 5.0):
            if c2 != c:
                res["eck"].append(zerlege(("S2", c), eck[("S2", c)], ("S2", c2)))
    for r in (0.003, 0.01, 0.03):
        for r2 in (0.003, 0.01, 0.03):
            res["fest"].append(zerlege(("S2", r), fest[("S2", r)], ("S1", r2)))
    for key, t in list(eck.items()) + list(fest.items()):
        j = tindex[key]
        gd = rho * horner(qt, t / rho) * horner(lam, rho / t)
        res["diag"].append({"t": str(key), "gen": [float(mp.re(gd)), float(mp.im(gd))], "nu": [float(mp.re(t * horner(nu, t / rho))), float(mp.im(t * horner(nu, t / rho)))],
                            "xi": [float(mp.re(horner(xis[j], t / rho))), float(mp.im(horner(xis[j], t / rho)))]})
    res["zeit"] = time.time() - t0
    return res


def c2(v):
    return complex(v[0], v[1])


def main():
    pfad = os.path.join(HIER, "k60.json")
    alle = json.load(open(pfad, encoding="utf-8")) if os.path.exists(pfad) else {}
    for a in sys.argv[1:]:
        r = lauf(int(a)); alle[a] = r
        json.dump(alle, open(pfad, "w", encoding="utf-8"), indent=1)
        print("n = %s fertig (%.0f s), max Prüfung (CD1) = %.1e" % (a, r["zeit"], max(p["pruef"] for p in r["eck"] + r["fest"])))
    ks = sorted(alle, key=int)
    print("\n(a) Eckskala: |XI|/|GEN| und |NU|/|GEN| (z auf Σ₂, t auf Σ₁ bzw. Σ₂; c = r/|d_n|) für n = " + " / ".join(ks))
    for i, p in enumerate(alle[ks[0]]["eck"]):
        zeile = []
        for k in ks:
            q = alle[k]["eck"][i]
            zeile.append("%.3f %.3f" % (abs(c2(q["xi"])) / abs(c2(q["gen"])), abs(c2(q["nu"])) / abs(c2(q["gen"]))))
        print("   z %-14s t %-14s  " % (p["z"], p["t"]) + "  |  ".join(zeile))
    print("\n(b) feste Punkte: XI, NU (komplex) je n")
    for i, p in enumerate(alle[ks[0]]["fest"]):
        print("   z %-16s t %-16s  XI: " % (p["z"], p["t"]) + "  ".join("%+.4f%+.4fi" % tuple(alle[k]["fest"][i]["xi"]) for k in ks)
              + "   NU: " + "  ".join("%+.4f%+.4fi" % tuple(alle[k]["fest"][i]["nu"]) for k in ks))
    print("\n(c) Diagonale t = z: GEN(t,t) = q_{n+1}(t)λ_n(t), NU, XI je n")
    for i, p in enumerate(alle[ks[0]]["diag"]):
        print("   t %-16s GEN: " % p["t"] + "  ".join("%+.4f%+.4fi" % tuple(alle[k]["diag"][i]["gen"]) for k in ks)
              + "   XI: " + "  ".join("%+.4f%+.4fi" % tuple(alle[k]["diag"][i]["xi"]) for k in ks))


if __name__ == "__main__":
    main()

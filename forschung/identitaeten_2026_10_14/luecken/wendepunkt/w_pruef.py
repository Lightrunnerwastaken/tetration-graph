"""Prüfungen der Eingabeannahme H (H1: Motorwerte auf 1e-30 genau; H2: |sexp| <= 10 auf |z - c| <= 0.3)
und des Kneser-Kennzeichens. Gleitkomma (mpmath), also Belege, kein Beweis.

P1  Genauigkeit: dieselben Knoten bei verschiedenen Motorgenauigkeiten (e: dps 130/100/70; 2: dps 80/60).
P2  Turmgleichung: exp(Motor auf |z-(c-1)|=0.2, dps 80) gegen Motor auf |z-c|=0.2 (dps 130), Knoten für Knoten (nur e).
P3  Kneser-Kennzeichen: sexp = F(z + θ(z)), θ = Σ_{n>=0} a_n e^{2πinz} (nur nicht-negative Moden, Koeffizienten aus
    etappe2, berechnet auf Im z = 0.003) gegen die Motorwerte an allen Kreisknoten mit Im z >= 0.1.
P4  Grobe Schranke H2: max |sexp| auf |z - c| = 0.3 (θ-Reihe für Im >= 0.1, sonst Taylor-Daten).
Aufruf: python w_pruef.py        (Basis 2 dps 60 einmal über den Motor, eigener Speicher luecken/wendepunkt/cache; ~1-3 min)
"""
import json
import sys

import mpmath as mp

from w_daten import E3, HIER, kreis

ZEN = ["-0.4", "-0.2", "0", "0.2", "0.4"]
mp.mp.dps = 140


def mpz(t):
    return mp.mpc(mp.mpf(t[0]), mp.mpf(t[1])), mp.mpc(mp.mpf(t[2]), mp.mpf(t[3]))


def p1():
    print("P1  Motorgenauigkeit an gleichen Knoten (max |Δ| je Kreis)")
    res = {}
    for c in ZEN:
        a = kreis("e", 130, c, "0.2", 256); b = kreis("e", 100, c, "0.2", 128); d = kreis("e", 70, c, "0.2", 256)
        d1 = max(abs(mpz(a[2 * j])[1] - mpz(b[j])[1]) for j in range(128))
        d2 = max(abs(mpz(a[j])[1] - mpz(d[j])[1]) for j in range(256))
        res[f"e {c}"] = [mp.nstr(d1, 3), mp.nstr(d2, 3)]
        print(f"   e  c={c:5s}  |dps130 - dps100| = {mp.nstr(d1, 3):9s}  |dps130 - dps70| = {mp.nstr(d2, 3)}")
    # Basis 2: dieselben Ausdrücke mit dps 60 über den Motor (eigener Speicher)
    sys.path.insert(0, str(E3.parent / "etappe2"))
    import kern
    kern.CDIR = HIER / "cache"
    for c in ZEN:
        a = kreis("2", 80, c, "0.2", 256)
        ausdr = [f"sexp(({t[0]}+({t[1]})*I))" if t[1] != "0" else f"sexp(({t[0]}))" for t in a]
        v = kern.werte("2", ausdr, 60)
        dd = max(abs(mpz(a[j])[1] - v[j]) for j in range(256))
        res[f"2 {c}"] = mp.nstr(dd, 3)
        print(f"   2  c={c:5s}  |dps80 - dps60| = {mp.nstr(dd, 3)}")
    return res


def p2():
    print("P2  Turmgleichung sexp(z+1) = exp(sexp z): Kreis um c-1 (dps 80) gegen Kreis um c (dps 130)")
    res = {}
    for c, cl in [("0.4", "-0.6"), ("0.2", "-0.8"), ("0", "-1.0"), ("-0.2", "-1.2"), ("-0.4", "-1.4")]:
        a = kreis("e", 130, c, "0.2", 256); b = kreis("e", 80, cl, "0.2", 256)
        diffs = []
        for j in range(256):
            za, sa = mpz(a[j]); zb, sb = mpz(b[j])
            assert abs(zb + 1 - za) < mp.mpf("1e-80")
            diffs.append(abs(mp.exp(sb) - sa))
        schlecht = [j for j in range(256) if diffs[j] > mp.mpf("1e-60")]
        gut = max(diffs[j] for j in range(256) if j not in schlecht)
        res[c] = dict(max_gut=mp.nstr(gut, 3), ausreisser=[(j, mp.nstr(diffs[j], 3), b[j][0][:12], b[j][1][:12]) for j in schlecht])
        print(f"   c={c:5s} max |Δ| = {mp.nstr(gut, 3)}   Ausreisser (> 1e-60): {res[c]['ausreisser']}")
    return res


class ThetaSexp:
    """sexp = Ψ(λ^{z+θ(z)}) für Basis b (λ = L ln b), Koeffizienten aus etappe2 (wie etappe3/s5_basis_b.py)."""

    def __init__(self, basis, K=70):
        self.b = mp.e if basis == "e" else mp.mpf(basis)
        lb = self.lb = mp.log(self.b)
        L = -mp.lambertw(-lb) / lb
        if mp.im(L) < 0:
            L = mp.conj(L)
        lam = self.lam = L * lb
        self.llam = mp.log(lam)
        P = [L, mp.mpc(1)] + [mp.mpc(0)] * (K - 1)
        for k in range(2, K + 1):
            q = [mp.mpc(0)] + [lb * x for x in P[1:k]] + [mp.mpc(0)] * (K + 1 - k)
            e = [mp.mpc(1)] + [mp.mpc(0)] * k
            for m in range(1, k + 1):
                e[m] = mp.fsum(j * q[j] * e[m - j] for j in range(1, m + 1)) / m
            P[k] = L * e[k] / (lam ** k - lam)
        self.P = P
        d = json.loads((E3.parent / "etappe2" / "daten" / f"a_{basis}_80_0.003_16384.json").read_text(encoding="utf-8"))
        self.A = [mp.mpc(mp.mpf(x[0]), mp.mpf(x[1])) for x in d["a"]]
        self.negativ = d.get("a_neg")

    def Psi(self, u):
        n = 0
        while abs(u) > mp.mpf("0.02"):
            u /= self.lam; n += 1
        w = mp.polyval(self.P[::-1], u)
        for _ in range(n):
            w = mp.exp(self.lb * w)
        return w

    def __call__(self, z):
        q = mp.expjpi(2 * z); th = mp.mpc(0); qn = mp.mpc(1)
        for ak in self.A:
            th += ak * qn; qn *= q
            if abs(qn) < mp.mpf(10) ** (-mp.mp.dps - 5):
                break
        return self.Psi(mp.exp(self.llam * (z + th)))


def p3_p4():
    res = {}
    for basis, dps in (("e", 130), ("2", 80)):
        mp.mp.dps = 60
        T = ThetaSexp(basis)
        print(f"P3  Basis {basis}: θ-Reihe (nur Moden n >= 0) gegen Motor, Knoten mit Im z >= 0.1")
        for c in ZEN:
            a = kreis(basis, dps, c, "0.2", 256)
            m = mp.mpf(0); n = 0
            for j in range(0, 256, 2):
                z, s = mpz(a[j])
                if mp.im(z) >= mp.mpf("0.1"):
                    m = max(m, abs(T(z) - s)); n += 1
            res[f"P3 {basis} {c}"] = mp.nstr(m, 3)
            print(f"   c={c:5s}  {n} Knoten, max |θ-sexp - Motor| = {mp.nstr(m, 3)}")
        # P4: |sexp| auf |z - c| = 0.3
        dat = json.loads((E3 / "daten" / ("taylor_130_256.json" if basis == "e" else "taylor_b2_80_256.json")).read_text(encoding="utf-8"))
        dat = dat.get("zentren", dat)
        co = {mp.mpf(k): [mp.mpf(x) for x in v] for k, v in dat.items()}
        mx = mp.mpf(0)
        for c in ZEN:
            cc = mp.mpf(c)
            for j in range(64):
                z = cc + mp.mpf("0.3") * mp.expjpi(mp.mpf(j) / 32)
                if abs(mp.im(z)) >= mp.mpf("0.1"):
                    v = T(z if mp.im(z) > 0 else mp.conj(z))
                else:
                    v = mp.polyval(co[cc][::-1], z - cc)
                mx = max(mx, abs(v))
        res[f"P4 {basis}"] = mp.nstr(mx, 6)
        print(f"P4  Basis {basis}: max |sexp| auf den Kreisen |z - c| = 0.3: {mp.nstr(mx, 6)}  (Annahme H2: <= 10)")
    return res


if __name__ == "__main__":
    out = {"P1": p1(), "P2": p2(), "P34": p3_p4()}
    (HIER / "pruefung.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")

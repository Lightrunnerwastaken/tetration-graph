"""Strang A, Runde 2: lineares Teilproblem P_N := Π_N f_K = J_N(f_K)⁻¹ T_N[eˣ + (f_K'∘f_K)·f_K] (Strang B, §4.1).

Exakt mit f'∘f = eˣ/f' (aus f∘f = exp):  J_N[k][j] = [f^j]_k + [eˣ/f']_{k−j},  rechte Seite T_N[eˣ + (eˣ/f')·f].
Benötigt nur a_0 … a_{N+1}. Koeffizientenquelle: DFT-Werte (gemeinsam/kneser_halbexp_koeff_1000.json) bis k ≤ KDFT,
darüber die Darboux-Reihe (gemeinsam/kneser_singularitaet_L.json, n ≤ 3, j ≤ 20; rel. 1e-21 … 1e-24 ab k = 500).
"""
import json
import sys
from pathlib import Path

import mpmath as mp
from flint import arb, arb_mat, arb_poly, ctx

HIER = Path(__file__).resolve().parent
GEM = HIER.parent / "gemeinsam"
DATEN = HIER / "daten"
KDFT = 700


def kneser_koeff(kmax, dps=70):
    """a_0..a_kmax als mpmath-Zahlen (DFT bis KDFT, danach Darboux); Zwischenspeicher daten/kneser_koeff_lang.json."""
    cache = DATEN / "kneser_koeff_lang.json"
    if cache.is_file():
        c = json.loads(cache.read_text(encoding="utf-8"))
        if len(c) > kmax:
            with mp.workdps(dps):
                return [mp.mpf(x) for x in c[:kmax + 1]]
    a = _kneser_koeff(kmax, dps)
    with mp.workdps(dps):
        cache.write_text(json.dumps([mp.nstr(x, 50) for x in a]), encoding="utf-8")
    return a


def _kneser_koeff(kmax, dps=70):
    with mp.workdps(dps):
        a = [mp.mpf(x) for x in json.loads((GEM / "kneser_halbexp_koeff_1000.json").read_text(encoding="utf-8"))["a"]]
        a = a[:min(KDFT, kmax) + 1]
        if kmax > KDFT:
            L = -mp.lambertw(-1)
            L = mp.conj(L) if mp.im(L) < 0 else L
            G = 2j * mp.pi / L
            LML = mp.log(-L)
            s = json.loads((GEM / "kneser_singularitaet_L.json").read_text(encoding="utf-8"))
            cx = lambda p: mp.mpc(mp.mpf(p[0]), mp.mpf(p[1]))  # noqa: E731
            terme = []
            for n in (1, 2, 3):
                dn = cx(s["d_n"][str(n)])
                e = [mp.mpc(1)] + [cx(p) for p in s["e_nj"][str(n)]]
                for j in range(21):
                    al = 1 + n * G + j
                    terme.append((dn * e[j] * mp.exp(al * LML), al))
            for k in range(KDFT + 1, kmax + 1):
                tot = mp.mpc(0)
                for c, al in terme:
                    tot += c * mp.exp(mp.loggamma(k - al) - mp.loggamma(-al) - mp.loggamma(k + 1))
                a.append(2 * mp.re(tot * L ** (-k)))
        return a


def pi_N(a_mp, N, prec, algorithm="approx", rauschen=None):
    """Koeffizienten b_0..b_N von Π_N f_K (arb, Mittelpunkte). a_mp: a_0..a_{N+1}."""
    alt = ctx.prec
    ctx.prec = prec
    try:
        a = [arb(mp.nstr(x, 75)) for x in a_mp[:N + 2]]
        if rauschen is not None:
            a = [x + arb(mp.nstr(r, 20)) for x, r in zip(a, rauschen)]
        f = arb_poly(a[:N + 1])
        fp = [(k + 1) * a[k + 1] for k in range(N + 1)]
        ex = [arb(1)]
        for k in range(1, N + 1):
            ex.append(ex[-1] / k)
        # E = eˣ/f' als Reihe bis x^N
        E = []
        for k in range(N + 1):
            s = ex[k]
            for i in range(1, k + 1):
                s -= fp[i] * E[k - i]
            E.append(s / fp[0])
        Ep = arb_poly(E)
        spalten = []
        P = arb_poly([1])
        for j in range(N + 1):
            c = P.coeffs()
            c = c + [arb(0)] * (N + 1 - len(c))
            spalten.append(c[:N + 1])
            P = (P * f).truncate(N + 1)
        zeilen = []
        for k in range(N + 1):
            zeilen.append([spalten[j][k] + (E[k - j] if j <= k else 0) for j in range(N + 1)])
        J = arb_mat(zeilen)
        rs = (Ep * f).truncate(N + 1)
        rc = rs.coeffs() + [arb(0)] * (N + 1)
        rhs = arb_mat([[ex[k] + rc[k]] for k in range(N + 1)])
        X = J.solve(rhs, algorithm=algorithm)
        return [arb(X[k, 0].mid()) for k in range(N + 1)]
    finally:
        ctx.prec = alt


def speichere(N, prec, b, tag=""):
    pfad = DATEN / f"lin_PN{tag}.json"
    d = json.loads(pfad.read_text(encoding="utf-8")) if pfad.is_file() else {}
    d[str(N)] = {"prec": prec, "b": [x.mid().str(60, radius=False) for x in b]}
    pfad.write_text(json.dumps(d), encoding="utf-8")


def lade(tag=""):
    pfad = DATEN / f"lin_PN{tag}.json"
    if not pfad.is_file():
        return {}
    return {int(k): v for k, v in json.loads(pfad.read_text(encoding="utf-8")).items()}


def newton_N(b_str, N, prec, schritte=2, stellen=60):
    """Newton für F_N(p) = T_N[p∘p − eˣ] = 0, Start p = Π_N f_K (Koeffizienten als Strings). Rückgabe: p (Strings),
    Liste der Korrekturgrössen max_k |d_k|·|L|^k."""
    alt = ctx.prec
    ctx.prec = prec
    try:
        with mp.workdps(60):
            L = -mp.lambertw(-1)
            RL = abs(L)
        a = [x if isinstance(x, arb) else arb(x) for x in b_str]
        korr = []
        for _ in range(schritte):
            f = arb_poly(a)
            spalten = []
            P = arb_poly([1])
            g, fpf = arb_poly([0]), arb_poly([0])
            for j in range(N + 1):
                c = P.coeffs()
                c = c + [arb(0)] * (N + 1 - len(c))
                spalten.append(c[:N + 1])
                g += a[j] * P
                if j + 1 <= N:
                    fpf += (j + 1) * a[j + 1] * P
                P = (P * f).truncate(N + 1)
            gc = g.coeffs() + [arb(0)] * (N + 1)
            fc = fpf.coeffs() + [arb(0)] * (N + 1)
            fak = arb(1)
            Fv = []
            for k in range(N + 1):
                fak = fak * k if k else fak
                Fv.append([gc[k] - 1 / fak])
            zeilen = [[spalten[j][k] + (fc[k - j] if j <= k else 0) for j in range(N + 1)] for k in range(N + 1)]
            d = arb_mat(zeilen).solve(arb_mat(Fv), algorithm="approx")
            a = [arb((a[i] - d[i, 0]).mid()) for i in range(N + 1)]
            with mp.workdps(30):
                korr.append(max(abs(mp.mpf(d[k, 0].mid().str(25, radius=False))) * RL ** k for k in range(N + 1)))
        return [x.mid().str(stellen, radius=False) for x in a], korr
    finally:
        ctx.prec = alt

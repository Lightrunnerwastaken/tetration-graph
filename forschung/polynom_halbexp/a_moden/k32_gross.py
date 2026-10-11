"""Strang A, Runde 2: sehr grosses N (Standard 1500) in Etappen (jede ≤ ~15 min).
  python k32_gross.py pi N       → Π_N f_K berechnen (lin_PN.json)
  python k32_gross.py newton N   → EIN Newton-Schritt vom letzten Iterat (daten/pN_gross_schritte.json; Start Π_N f_K)
  python k32_gross.py bericht N  → Kenngrössen: |c_N|, K_N, |P(L) − L|, max|P| auf D, Newton-Korrekturen, |p − P|
"""
import json
import sys
import time

sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402

from lin_kern import DATEN, kneser_koeff, lade, newton_N, pi_N, speichere  # noqa: E402

mp.mp.dps = 60
modus, N = sys.argv[1], int(sys.argv[2])
prec = int(0.72 * N) + 200
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
RL = abs(L)
pfad = DATEN / "pN_gross_schritte.json"
t0 = time.time()
if modus == "pi":
    a = kneser_koeff(N + 1)
    b = pi_N(a, N, prec)
    speichere(N, prec, b)
    print(f"Π_{N} f_K fertig ({time.time() - t0:.0f} s, prec {prec})")
elif modus == "newton":
    st = json.loads(pfad.read_text(encoding="utf-8")) if pfad.is_file() else {}
    eintrag = st.get(str(N), {"iterate": [], "korr": []})
    start = eintrag["iterate"][-1] if eintrag["iterate"] else lade()[N]["b"]
    p, korr = newton_N(start, N, prec, schritte=1)
    eintrag["iterate"].append(p)
    eintrag["korr"].append(mp.nstr(korr[0], 4))
    st[str(N)] = eintrag
    pfad.write_text(json.dumps(st), encoding="utf-8")
    print(f"N = {N}: Newton-Schritt {len(eintrag['korr'])}, Korrektur {mp.nstr(korr[0], 4)} ({time.time() - t0:.0f} s)")
else:
    import numpy as np
    b = [mp.mpf(x) for x in lade()[N]["b"]]
    a = kneser_koeff(4 * N)
    ld = json.loads((DATEN / "linie48.json").read_text(encoding="utf-8"))
    ws = [mp.mpf(x) for x in ld["w"]]
    zs = [mp.mpf(x) for x in ld["z"]]
    Tp = [mp.mpf(x) for x in ld["Tp"]]
    hv = lambda c, z: mp.fsum(x * z ** k for k, x in enumerate(c))  # noqa: E731
    fz = [hv(a[:1001], z) for z in zs]
    th_ = [(hv(b, z) - f) / t for z, f, t in zip(zs, fz, Tp)]
    c = mp.fsum(t * mp.expjpi(2 * w) for t, w in zip(th_, ws)) / len(ws)
    G = 2j * mp.pi / L
    sing = json.loads((DATEN.parent.parent / "gemeinsam" / "kneser_singularitaet_L.json").read_text(encoding="utf-8"))
    cx = lambda p: mp.mpc(mp.mpf(p[0]), mp.mpf(p[1]))  # noqa: E731
    D1 = cx(sing["d_n"]["1"]) * mp.exp((1 + G) * mp.log(-L))
    k0 = cx(sing["kappa"]["0"])
    Dd = L * mp.exp(L / 2) * mp.exp(-2j * mp.pi * k0) * mp.exp((1 - G) * mp.log(-L))
    koef = lambda al: mp.exp(mp.loggamma(N - al) - mp.loggamma(-al) - mp.loggamma(N + 1)) * L ** (-N)  # noqa: E731
    K = c * Dd * koef(1 - G) / (N * D1 * koef(1 + G))
    PL = hv(b, L) - L
    cc = np.array([float(x * RL ** k) for k, x in enumerate(b)])
    Mg = 1 << int(np.ceil(np.log2(16 * (N + 1))))
    mx = float(np.max(np.abs(np.fft.ifft(cc, Mg) * Mg)))
    print(f"N = {N}: |c_N| = {mp.nstr(abs(c), 5)}, arg {mp.nstr(mp.arg(c), 5)}, K_N = {mp.nstr(K, 6)}; |P(L) − L| = {mp.nstr(abs(PL), 5)}; "
          f"max_D|P_N| = {mx:.6f}")
    if pfad.is_file() and str(N) in json.loads(pfad.read_text(encoding="utf-8")):
        e = json.loads(pfad.read_text(encoding="utf-8"))[str(N)]
        p = [mp.mpf(x) for x in e["iterate"][-1]]
        d = [x - y for x, y in zip(p, b)]
        print(f"   Newton-Korrekturen {e['korr']}; |p_N − P_N|(L) = {mp.nstr(abs(hv(d, L)), 4)}, |E(L)|² = {mp.nstr(abs(PL) ** 2, 4)}; "
              f"p_N(1) − T(½) = {mp.nstr(mp.fsum(p) - mp.fsum(a[:4 * N + 1]), 5)}")

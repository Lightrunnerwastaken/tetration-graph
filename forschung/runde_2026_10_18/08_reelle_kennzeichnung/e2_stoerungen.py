"""E2: gestörte Tetrationen T = sexp_K(x + p(x)): erste negative ungerade Ableitung, Wurzeltest-Radius.
Aufruf: python e2_stoerungen.py N dps name1 name2 ...   (Gitter x = -1.5 … 0.4, Schritt 0.1)"""
import sys, json, mpmath as mp, rk_lib as R, stoerungen as S
N = int(sys.argv[1]); dps = int(sys.argv[2]); namen = sys.argv[3:]
R.init(dps)
out = {}
for name in namen:
    st = S.katalog(name)
    zeilen = []
    for i in range(20):
        x = mp.mpf('-1.5') + mp.mpf(i) / 10
        try:
            c = R.tet(st, x, N)
        except ValueError as e:
            continue
        neg = [k for k in range(1, N, 2) if c[k] <= 0]
        # Wurzeltest-Radius aus k in [N-20, N)
        rs = [abs(c[k]) ** (-mp.mpf(1) / k) for k in range(N - 20, N) if c[k] != 0]
        rho = min(rs)
        zeilen.append({"x": float(x), "erste_neg_ungerade": neg[0] if neg else None, "anzahl_neg_ungerade": len(neg),
                       "rho_wurzeltest": float(rho), "x_plus_2": float(x + 2)})
        print(name, mp.nstr(x, 3), "erste neg. ungerade:", neg[0] if neg else "-", "(", len(neg), ")",
              " rho~", mp.nstr(rho, 4), " x+2 =", mp.nstr(x + 2, 3), flush=True)
    out[name] = zeilen
json.dump(out, open(f"daten/e2_{'_'.join(n.replace(':','') for n in namen)}_N{N}.json", "w"), indent=1)

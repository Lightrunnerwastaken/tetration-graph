# t-Koordinate: Basis b mit ln b = x = t e^{-t} (t = ln h(b), h = unendlicher Turm).
# tau_n(t) = eps_n(t e^{-t}) = ln(b^^(n+1)),  tau_n = t + delta_n,  delta_{k+1} = t (e^{delta_k} - 1),  delta_0 = t(e^{-t}-1).
# (A) Untere Schranke L_n = max_r (r - M_n(r)), M_n(r) = max_{|t|=r} |tau_n(t) - t|  (Rouche => R_n >= L_n).
# (B) Lage der sichtbaren Verzweigungspunkte in t: |t_c| fuer die in radius.py gefundenen c.
# (C) Reelle kritische Punkte (n ungerade) = Minima der geraden Tuerme y^^(n+1) auf (0,1).
import sys, cmath, math
import mpmath as mp

def tau(n, t):
    d = t * (cmath.exp(-t) - 1)
    for _ in range(n):
        d = t * (cmath.exp(d) - 1)
    return t + d

def M(n, r, m=720):
    return max(abs(tau(n, r * cmath.exp(2j * math.pi * k / m)) - r * cmath.exp(2j * math.pi * k / m)) for k in range(m))

def L(n):
    best = (-1, None)
    for i in range(1, 200):
        r = i / 200
        try:
            val = r - M(n, r)
        except OverflowError:
            break
        if val > best[0]:
            best = (val, r)
    return best

R = {1: 0.367879441171442, 2: 0.449544113030792, 3: 0.517584939137185, 4: 0.559531199099174, 5: 0.600241098018159,
     6: 0.630976549641869, 7: 0.655752155914705, 8: 0.68089215989053, 10: 0.714867528693796, 12: 0.741978208982033,
     16: 0.783215381779534}
C = {2: -0.0372507656791 + 0.819937714011j, 4: -0.608201181522 + 0.979510361332j, 6: -0.985996202264 + 0.957640102245j,
     8: 0.12409141892 + 1.11074360429j, 10: -0.185576702657 + 1.29359290956j, 12: -0.465350111017 + 1.37859687105j,
     16: -0.918825676109 + 1.40200874347j}

if __name__ == "__main__":
    teil = sys.argv[1] if len(sys.argv) > 1 else "ABC"
    if "A" in teil:
        print("(A) n, R_n, L_n (untere Schranke), r*, n(1-R_n)/ln n")
        for n in sorted(R):
            l, r = L(n)
            print(f"  n={n:3d}  R={R[n]:.6f}  L={l:.6f}  r*={r}  n(1-R)/ln n={(n*(1-R[n])/math.log(n)) if n>1 else float('nan'):.4f}")
    if "B" in teil:
        print("(B) sichtbarer Verzweigungspunkt c in t-Koordinate: t = -W(-c)")
        for n, c in sorted(C.items()):
            t = complex(-mp.lambertw(-c))
            print(f"  n={n:3d}  c={c:.6f}  t={t:.6f}  |t|={abs(t):.6f}  arg t/2pi={cmath.phase(t)/(2*math.pi):.5f}  tau_n(t)={tau(n,t):.6f}")
    if "C" in teil:
        mp.mp.dps = 30
        print("(C) reelle kritische Punkte c_n (n ungerade): Minimum von y^^(n+1) bei y=e^{c_n}, Wert exp(eps_n(c_n))")
        def ed2(n, x):
            e, d, dd = x, mp.mpf(1), mp.mpf(0)
            for _ in range(n):
                E = mp.exp(e)
                e, d, dd = x * E, E + x * E * d, 2 * E * d + x * E * (d * d + dd)
            return e, d, dd
        z = mp.mpf(-1)
        for n in list(range(1, 40, 2)) + [51, 81, 101, 151, 201, 301, 401]:
            # Startwert: vorheriger Wert, leicht nach links
            for start in [z, z - mp.mpf('0.05'), z - mp.mpf('0.2')]:
                w = start
                try:
                    for _ in range(100):
                        e, d, dd = ed2(n, w)
                        st = d / dd
                        w -= st
                        if abs(st) < mp.mpf(10) ** -22:
                            break
                    if abs(ed2(n, w)[1]) < 1e-15 and -mp.e < w < z + mp.mpf('1e-9'):
                        break
                except Exception:
                    continue
            z = w
            v = ed2(n, z)[0]
            print(f"  n={n:4d}  c={mp.nstr(z,15)}  y=e^c={mp.nstr(mp.exp(z),12)}  min y^^(n+1)={mp.nstr(mp.exp(v),12)}  (c+e)*sqrt(n)={mp.nstr((z+mp.e)*mp.sqrt(n),8)}  (v+1)*sqrt(n)={mp.nstr((v+1)*mp.sqrt(n),8)}")

"""E3: komplexe Lösungen von φ(z) = −2 (Log-Verzweigung von sexp_K∘φ) in der Halbebene Re z > −2, und
ob das Bild der Strecke [x, z0] unter φ den Schnitt (−∞, −2] vorher kreuzt (dann anderes Blatt, Singularität unsicher).
Aufruf: python e3_singularitaeten.py name ..."""
import sys, mpmath as mp, stoerungen as S
mp.mp.dps = 30
for name in sys.argv[1:]:
    st = S.katalog(name)
    f = lambda z: z + S.p_komplex(st, z)[0] + 2
    df = lambda z: 1 + S.p_komplex(st, z)[1]
    wurzeln = []
    for xr in [i / 4 for i in range(-8, 9)]:
        for yi in [0.2 * j for j in range(1, 16)]:
            try:
                z = mp.findroot(f, mp.mpc(xr, yi), df=df, tol=1e-25, maxsteps=60)
            except Exception:
                continue
            if mp.im(z) > 1e-8 and mp.re(z) > -2 and abs(f(z)) < 1e-20 and all(abs(z - w) > 1e-8 for w in wurzeln):
                wurzeln.append(z)
    wurzeln.sort(key=lambda z: mp.re(z))
    print(name, "Wurzeln von φ = −2 mit Re > −2, 0 < Im ≤ 3:")
    for z in wurzeln:
        # Blattprüfung: Strecke von Re z (reell) senkrecht nach z
        x0 = mp.re(z); kreuzt = False
        for t in range(1, 400):
            w = mp.mpc(x0, mp.im(z) * t / 400)
            v = w + S.p_komplex(st, w)[0]
            v2 = mp.mpc(x0, mp.im(z) * (t + 1) / 400); v2 = v2 + S.p_komplex(st, v2)[0]
            if t < 399 and (mp.im(v) > 0) != (mp.im(v2) > 0) and min(mp.re(v), mp.re(v2)) < -2:
                kreuzt = True
        print("   z0 =", mp.nstr(z, 8), " |Im| =", mp.nstr(mp.im(z), 5), "  senkrechter Weg kreuzt Schnitt:", kreuzt)

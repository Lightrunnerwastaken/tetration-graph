"""Sucht komplexe Fixpunkte von sexp_2, die fuer slog (Hauptzweig) anziehend sind: alle 9 Startpunkte laufen nach p = -1.7439."""
import mpmath as mp, motor
mp.mp.dps=30
zs=[mp.mpc(*s) for s in [(0.5,0.5),(-1.9,1.0),(-1.5,2),(0.5,3),(1,1.5),(-1.95,0.3),(2,0.5),(-1.99,1.4),(1.47,1.66)]]
for it in range(150):
    zs=[motor.roh("2",[f"iferr(slog({motor.gpzahl(z,35)}),E,1000)"],30)[0] if abs(z)<100 else z for z in zs]
for z in zs:
    d = motor.roh("2",[f"iferr(sexp({motor.gpzahl(z,35)}),E,1000)"],30)[0]
    print(mp.nstr(z,15), " |sexp(z)-z|=", mp.nstr(abs(d-z),3))

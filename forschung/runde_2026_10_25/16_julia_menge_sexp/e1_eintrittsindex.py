"""Eintrittsindex N(z) = min{n : Re F(z+n) > RG} (F = Kneser-sexp Basis e) und Streifenaustritt T(z) = min{n : |Im F(z+n)| > pi}.
Beides gilt exakt N(z+1) = N(z) - 1 (Verschiebung), daher reicht der Grundstreifen 0 <= Re z < 1.
Prueft: Steigung von N in Im z gegen Im L / Re L, Gitterperiode omega = 2 pi i / L - 4, Lage von F(Streifen) zu S_0.
Bilder: bilder/e1_zebene.png, bilder/e1_webene.png.  Ausgabe: cache/e1_zahlen.json
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from kneser_feld import feld

HIER = Path(__file__).resolve().parent
(HIER / "bilder").mkdir(exist_ok=True)
L = complex(0.31813150520476413, 1.3372357014306895)
RG = 50.0
NMAX = 90


def indizes(w):
    """Vorwaerts-Bahn von exp: erster Eintritt Re>RG und erster Austritt |Im|>pi (NMAX = nie)."""
    w = w.copy()
    N = np.full(w.shape, NMAX, dtype=int)
    T = np.full(w.shape, NMAX, dtype=int)
    aktiv = np.ones(w.shape, bool)
    for n in range(NMAX):
        neuN = aktiv & (w.real > RG)
        N[neuN] = n
        neuT = (T == NMAX) & (np.abs(w.imag) > np.pi) & aktiv
        T[neuT] = n
        aktiv &= ~neuN
        if not aktiv.any():
            break
        with np.errstate(over="ignore", invalid="ignore"):
            w[aktiv] = np.exp(w[aktiv])
    return N, T


def grundstreifen(nx, ny, ymax):
    x = (np.arange(nx) + 0.5) / nx
    y = np.linspace(-ymax, ymax, ny)
    Z = x[None, :] + 1j * y[:, None]
    return x, y, Z, feld(Z)


if __name__ == "__main__":
    out = {}
    # --- Grundstreifen fein
    x, y, Z, FZ = grundstreifen(200, 1201, 7.5)
    N, T = indizes(FZ)
    out["max_abs_Im_F_streifen"] = float(np.abs(FZ.imag).max())
    out["F_streifen_in_S0"] = bool(np.abs(FZ.imag).max() < np.pi)
    # --- Steigung: Median von N je Zeile gegen Im z (oben, y>1)
    oben = (y > 1.0) & (y < 7.4)
    med = np.array([np.median(N[i][N[i] < NMAX]) if (N[i] < NMAX).any() else np.nan for i in range(len(y))])
    mittel = np.array([np.mean(N[i][N[i] < NMAX] + x[N[i] < NMAX]) for i in range(len(y))])
    p = np.polyfit(y[oben], mittel[oben], 1)
    out["steigung_N_pro_Im"] = float(p[0])
    out["ImL_durch_ReL"] = L.imag / L.real
    out["anteil_nie_eingetreten"] = float((N == NMAX).mean())
    # --- Gitterperiode omega: N(z+omega) == N(z)+4 ?
    P = 2j * np.pi / L
    om = P - 4
    out["P"] = [P.real, P.imag]
    out["omega"] = [om.real, om.imag]
    zs = Z + om
    k = np.floor(zs.real)
    z0 = zs - k
    N2, _ = indizes(feld(z0))
    N2 = N2 - k.astype(int)
    gleich = (N2 == N + 4) & (N < NMAX - 10)
    gueltig = (N < NMAX - 10) & (N2 < NMAX)
    baender = []
    for lo in [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]:
        m = (y >= lo) & (y < lo + 0.5)
        g = gueltig[m]
        baender.append([lo, float(gleich[m][g].mean()) if g.any() else None])
    out["anteil_N(z+omega)=N(z)+4_je_Imband"] = baender
    # Gegenprobe mit anderer Schwelle: gleiche Statistik fuer T
    Tz = T
    T2 = indizes(feld(z0))[1] - k.astype(int)
    gT = (Tz < NMAX - 10) & (T2 < NMAX)
    bT = []
    for lo in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]:
        m = (y >= lo) & (y < lo + 0.5)
        bT.append([lo, float((T2 == Tz + 4)[m][gT[m]].mean())])
    out["anteil_T(z+omega)=T(z)+4_je_Imband"] = bT
    (HIER / "cache" / "e1_zahlen.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))

    # --- Bild z-Ebene: Re z in [-1,3] ueber Verschiebung, Groesse  N(z) + Re(Lz)/Re L  (asympt. Gitter-periodisch)
    xs = np.linspace(-1, 3, 800, endpoint=False)
    kk = np.floor(xs).astype(int)
    col = ((xs - kk) * 200).astype(int)
    Nbig = N[:, col] - kk[None, :]
    Tbig = T[:, col] - kk[None, :]
    Zb = xs[None, :] + 1j * y[:, None]
    fig, ax = plt.subplots(1, 2, figsize=(13, 9), sharey=True)
    A = np.where(N[:, col] < NMAX, Nbig % 6, np.nan)
    ax[0].imshow(A, origin="lower", extent=[-1, 3, -7.5, 7.5], cmap="viridis", aspect="auto", interpolation="nearest")
    ax[0].set_title("Eintrittsindex N(z) mod 6, N(z)=min{n: Re F(z+n)>50}\n(F = Kneser-sexp, Basis e; schwarz: kein Eintritt bis n=90)")
    B = np.where(T[:, col] < NMAX, Tbig % 6, np.nan)
    ax[1].imshow(B, origin="lower", extent=[-1, 3, -7.5, 7.5], cmap="magma", aspect="auto", interpolation="nearest")
    ax[1].set_title("Streifenaustritt T(z) mod 6, T(z)=min{n: |Im F(z+n)|>pi}")
    for a in ax:
        a.set_facecolor("black")
        a.set_xlabel("Re z")
        for j in range(-1, 4):
            a.axvline(j, color="white", lw=0.4, alpha=0.5)
        # Gittervektoren
        z0p = 0.0 + 3.0j
        for m_ in range(4):
            q = z0p + m_ * om
            a.plot([q.real], [q.imag], "w+", ms=8)
    ax[0].set_ylabel("Im z")
    fig.tight_layout()
    fig.savefig(HIER / "bilder" / "e1_zebene.png", dpi=110)

    # --- Bild w-Ebene: Eintrittsindex von exp selbst + Bilder senkrechter Linien unter F
    u = np.linspace(-3, 5, 900)
    v = np.linspace(-4.5, 4.5, 900)
    W = u[None, :] + 1j * v[:, None]
    Nw, Tw = indizes(W)
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.imshow(np.where(Nw < NMAX, Nw % 6, np.nan), origin="lower", extent=[-3, 5, -4.5, 4.5], cmap="viridis", alpha=0.85)
    ax.set_facecolor("black")
    yy = np.linspace(-7.5, 7.5, 3001)
    for x0, c in [(0.0, "red"), (0.25, "orange"), (0.5, "yellow"), (0.75, "cyan"), (1.0, "red")]:
        cz = feld(x0 + 1j * yy) if x0 < 1 else np.exp(feld(0.0 + 1j * yy))
        ax.plot(cz.real, cz.imag, color=c, lw=1.4, label=f"F({x0}+iy)")
    for y0 in [0.25, 0.5, 1.0, 2.0]:
        xx = np.linspace(0, 1, 400, endpoint=False)
        cz = feld(xx + 1j * y0)
        ax.plot(cz.real, cz.imag, color="white", lw=0.8, ls="--")
        ax.plot(cz.real, -cz.imag, color="white", lw=0.8, ls="--")
    ax.plot([L.real, L.real], [L.imag, -L.imag], "w*", ms=12)
    ax.axhline(np.pi, color="white", lw=0.5)
    ax.axhline(-np.pi, color="white", lw=0.5)
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title("Dynamische Ebene von exp: Eintrittsindex mod 6 (Re>50),\nBilder senkrechter Linien unter Kneser-F (rot: Re z=0 und 1), weiss gestrichelt: Im z = 0.25, 0.5, 1, 2")
    fig.tight_layout()
    fig.savefig(HIER / "bilder" / "e1_webene.png", dpi=110)

"""Strang V — Ks E4′ (k_kern/BERICHT §28.5, PLAN „Anpassung K 10:48“): Q_n = max_{z ∈ K_{r₀}, t ∈ 𝒞} |ΔK_n(z, t)| für alle n ≤ N.

ΔK_n(z, t) = δ_n(z)·Λ_n(t) (Rang 1, Satz K5): δ_n(z) = ε_zᵀU_n⁻¹e_n, Λ_n(t) = e_nᵀL_n⁻¹ẽ_t aus EINER LU von J̃_N (wie v48), für viele Punkte zugleich:
  Zeilenblock r = ε_zᵀU⁻¹ (Rekursion r ← r − w_k·(Zeile k von U), Rang-1 in C), Spaltenblock b = L⁻¹ẽ_t (Vorwärts-Substitution, Rang-1 in C).
Mengen (obere Hälfte; die untere ist konjugiert, J̃ reell):
  Z: ∂K_{r₀} für r₀ = 0.1/0.05/0.02 (Kreis |z| = ρ ohne U_{r₀}(L), dazu Bogen |z − L| = r₀ in D̄), z = L (Normierung), ∂D-Punkte ρe^{i·120°/150°/180°}.
  T (Kontur 𝒞): Σ₁ (u_A = 60°, 80°) und Σ₂ (u_A = −90°): θ(r) = u_A + 71.46° + 1.623°·ln(r/ρ), r ∈ [1e-7, r₁], r₁ = 0.05 (Ks §28.2);
     fern: Γout-Bogen |t| = ρ mit |t − L| ≥ r₁, Linse |t| = ρ(1 + 0.03) über Γin mit |t − L| ≥ r₁; Verbindungsbögen |t − L| = r₁ vom Spiralende zum fernen Teil.
Ausgabe je n: |P_n|, R_P(r₀) = max_{∂K_{r₀}}|δ_n/δ_n(L)|, R_A je Teil (mit Lage), Q_n(r₀) = |P_n|·R_P(r₀)·R_A, ln|δ_n/δ_n(L)| bei 120°/150°/180°.
Aufruf: python v51_E4.py kubisch|kneser N [P]
"""
import cmath
import json
import math
import sys
import time
from pathlib import Path

HIER = Path(__file__).resolve().parent
MODELL, N = sys.argv[1], int(sys.argv[2])
P = int(sys.argv[3]) if len(sys.argv) > 3 else int(0.5 * N) + 300
sys.argv = ["v48", MODELL, str(N), str(P)]
sys.path.insert(0, str(HIER))
import v48_RL as V48  # noqa: E402
from flint import arb, acb, arb_mat, arb_poly, ctx  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
Lc = complex(float(V48.L_RE), float(V48.L_IM))
rho, th = abs(Lc), cmath.phase(Lc)
R1 = 0.05
TS = 1.623 * math.pi / 180                     # Spiralsteigung je e-Faktor


def spirale(uA, npkt=30):
    out = []
    for i in range(npkt):
        r = 1e-7 * (R1 / 1e-7) ** (i / (npkt - 1))
        ang = math.radians(uA + 71.46) + TS * math.log(r / rho)
        out.append(Lc + r * cmath.exp(1j * ang))
    return out


def punkte():
    Z, zname = [Lc], ["L"]
    for k in range(240):                         # Kreis, obere Hälfte (einmal; je r₀ nach Abstand gefiltert)
        phi = math.pi * k / 239
        z = rho * cmath.exp(1j * phi)
        if abs(z - Lc) >= 0.02:
            Z.append(z)
            zname.append(f"kreis:{math.degrees(phi):.1f}:{abs(z - Lc):.5f}")
    for r0 in (0.1, 0.05, 0.02):
        for k in range(40):                      # Bogen |z − L| = r₀ in D̄
            a = th + math.pi / 2 + math.pi * k / 39
            z = Lc + r0 * cmath.exp(1j * a)
            if abs(z) <= rho * (1 + 1e-12):
                Z.append(z)
                zname.append(f"K{r0}:bogen:{math.degrees(a) % 360:.1f}")
    for deg in (120, 150, 180):
        Z.append(rho * cmath.exp(1j * math.radians(deg)))
        zname.append(f"rand{deg}")
    T, tname = [Lc], ["L"]
    for uA in (60, 80, -90):
        for t in spirale(uA):
            T.append(t)
            tname.append(f"sigma{uA}")
    for k in range(100):                         # Γout-Bogen
        phi = th * k / 99
        t = rho * cmath.exp(1j * phi)
        if abs(t - Lc) >= R1:
            T.append(t)
            tname.append(f"fern:gout:{math.degrees(phi):.1f}")
    for k in range(100):                         # Linse s = 0.03 über Γin
        phi = th + (math.pi - th) * k / 99
        t = rho * 1.03 * cmath.exp(1j * phi)
        if abs(t - Lc) >= R1:
            T.append(t)
            tname.append(f"fern:linse:{math.degrees(phi):.1f}")
    # Verbindungsbögen |t − L| = r₁: von Σ₁(80°)-Ende zur Linse, von Σ₂-Ende zum Γout-Bogen
    a1 = math.radians(80 + 71.46) + TS * math.log(R1 / rho)
    a2 = math.radians(-90 + 71.46) + TS * math.log(R1 / rho)
    # Richtung des Linsen- bzw. Γout-Punkts im Abstand r₁ von L
    def richtung(rad):
        best = None
        for k in range(20001):
            phi = th - 0.2 + 0.4 * k / 20000
            t = rad * cmath.exp(1j * phi)
            d = abs(abs(t - Lc) - R1)
            if (rad > rho and phi > th) or (rad == rho and phi < th):
                if best is None or d < best[0]:
                    best = (d, cmath.phase(t - Lc))
        return best[1]
    al = richtung(rho * 1.03)
    ag = richtung(rho)
    while ag - a2 > math.pi:
        ag -= 2 * math.pi
    while a2 - ag > math.pi:
        ag += 2 * math.pi
    while al - a1 > math.pi:
        al -= 2 * math.pi
    while a1 - al > math.pi:
        al += 2 * math.pi
    print(f"  Verbindungsbögen: Linse {math.degrees(a1):.1f}° → {math.degrees(al):.1f}°, Γout {math.degrees(a2):.1f}° → {math.degrees(ag):.1f}°", flush=True)
    for k in range(12):
        T.append(Lc + R1 * cmath.exp(1j * (a1 + (al - a1) * k / 11)))
        tname.append("verb:linse")
    for k in range(12):
        T.append(Lc + R1 * cmath.exp(1j * (a2 + (ag - a2) * k / 11)))
        tname.append("verb:gout")
    return Z, zname, T, tname


def main():
    t0 = time.time()
    Z, zname, T, tname = punkte()
    ctx.prec = P
    fk, mk = V48.modell(N)
    f = arb_poly(fk)
    rh = (arb(V48.L_RE) ** 2 + arb(V48.L_IM) ** 2).sqrt()
    rp = [rh ** k for k in range(N + 1)]
    A = arb_mat(N + 1, N + 1)
    col = arb_poly([1])
    for j in range(N + 1):
        for k in range(N + 1):
            ck = col[k] if k <= col.degree() else arb(0)
            mm = mk[k - j] if k >= j else arb(0)
            A[k, j] = ((ck + mm) * rp[k] / rp[j]).mid()
        col = (col * f).truncate(N + 1)
    nz, nt = len(Z), len(T)
    print(f"[{MODELL}, N = {N}, P = {P} Bit] J̃ aufgebaut, {nz} z-Punkte, {nt} t-Punkte ({time.time() - t0:.0f}s)", flush=True)
    # r = ε_zᵀ: Zeilen 2i (Re), 2i+1 (Im); b = ẽ_t: Spalten 2j (Re), 2j+1 (Im). Potenzen in acb mit voller Präzision.
    r = arb_mat(2 * nz, N + 1)
    for i, z in enumerate(Z):
        q = acb(arb(z.real), arb(z.imag)) / acb(rh)
        if i == 0:
            q = acb(arb(V48.L_RE), arb(V48.L_IM)) / acb(rh)
        qk = acb(1)
        for k in range(N + 1):
            r[2 * i, k] = qk.real.mid()
            r[2 * i + 1, k] = qk.imag.mid()
            qk = qk * q
    b = arb_mat(N + 1, 2 * nt)
    for j, t in enumerate(T):
        q = acb(rh) / acb(arb(t.real), arb(t.imag))
        if j == 0:
            q = acb(rh) / acb(arb(V48.L_RE), arb(V48.L_IM))
        qk = acb(1)
        for k in range(N + 1):
            b[k, 2 * j] = qk.real.mid()
            b[k, 2 * j + 1] = qk.imag.mid()
            qk = qk * q
    print(f"  Startvektoren ({time.time() - t0:.0f}s)", flush=True)
    off, m = 0, N + 1
    erg = []
    for k in range(N + 1):
        i = k - off
        p = A[i, i]
        # w_k = r[:, i]/p (Spaltenvektor), Kopf δ_k(z) = w_k
        wv = arb_mat(2 * nz, 1)
        for ii in range(2 * nz):
            wv[ii, 0] = (r[ii, i] / p).mid()
        # y_k = b[i, :] (Zeile), Kopf Λ_k(t) = y_k
        yrow = arb_mat(1, 2 * nt)
        for jj in range(2 * nt):
            yrow[0, jj] = b[i, jj]
        # Auswertung für n = k
        dz = [complex(float(wv[2 * a, 0].mid()), float(wv[2 * a + 1, 0].mid())) for a in range(nz)]
        lt = [complex(float(yrow[0, 2 * c].mid()), float(yrow[0, 2 * c + 1].mid())) for c in range(nt)]
        if k >= 50:
            dL, lL = dz[0], lt[0]
            Pn = dL * lL
            rec = {"n": k, "P": abs(Pn)}
            for r0 in (0.1, 0.05, 0.02):
                vals = [(abs(dz[a]) / abs(dL), zname[a]) for a in range(nz)
                        if zname[a].startswith(f"K{r0}:") or (zname[a].startswith("kreis:") and float(zname[a].split(":")[2]) >= r0)]
                mx = max(vals)
                rec[f"RP{r0}"] = mx[0]
                rec[f"RP{r0}_lage"] = mx[1]
            for deg in (120, 150, 180):
                a = zname.index(f"rand{deg}")
                rec[f"lnrand{deg}"] = math.log(abs(dz[a]) / abs(dL)) if dz[a] != 0 else -1e9
            teile = {"sigma": [], "fern": [], "verb": []}
            for c in range(1, nt):
                key = tname[c].split(":")[0]
                key = "sigma" if key.startswith("sigma") else key
                teile[key].append((abs(lt[c]) / abs(lL), tname[c]))
            for key, vals in teile.items():
                mx = max(vals)
                rec[f"RA_{key}"] = mx[0]
                rec[f"RA_{key}_lage"] = mx[1]
            RA = max(rec["RA_sigma"], rec["RA_fern"], rec["RA_verb"])
            for r0 in (0.1, 0.05, 0.02):
                rec[f"Q{r0}"] = abs(Pn) * rec[f"RP{r0}"] * RA
            erg.append(rec)
            if k % 200 == 0:
                print(f"n = {k}: |P_n| = {abs(Pn):.5f}; R_P(0.1/0.05/0.02) = {rec['RP0.1']:.3g}/{rec['RP0.05']:.3g}/{rec['RP0.02']:.3g} "
                      f"({rec['RP0.02_lage']}); R_A Σ/fern/Verb = {rec['RA_sigma']:.4f}/{rec['RA_fern']:.3g} ({rec['RA_fern_lage']})/{rec['RA_verb']:.3g}; "
                      f"Q_n(0.05) = {rec['Q0.05']:.4f}, Q/|P| = {rec['Q0.05'] / abs(Pn):.4f}; ln|δ/δ(L)| 120/150/180° = "
                      f"{rec['lnrand120']:.2f}/{rec['lnrand150']:.2f}/{rec['lnrand180']:.2f}  ({time.time() - t0:.0f}s)", flush=True)
        # r ← r − w ⊗ (Zeile i von U), A und b Rang-1
        cvec = arb_mat(m, 1)
        rvec = arb_mat(1, m)
        for ii in range(i + 1, m):
            cvec[ii, 0] = (A[ii, i] / p).mid()
        for jj in range(i + 1, m):
            rvec[0, jj] = A[i, jj]
        r = (r - wv * rvec).mid()
        A = (A - cvec * rvec).mid()
        b = (b - cvec * yrow).mid()
        if (k + 1) % 250 == 0 and k < N:
            noff = k + 1
            mn = N + 1 - noff
            d = noff - off
            A2, b2, r2 = arb_mat(mn, mn), arb_mat(mn, 2 * nt), arb_mat(2 * nz, mn)
            for ii in range(mn):
                for jj in range(mn):
                    A2[ii, jj] = A[ii + d, jj + d]
                for jj in range(2 * nt):
                    b2[ii, jj] = b[ii + d, jj]
            for ii in range(2 * nz):
                for jj in range(mn):
                    r2[ii, jj] = r[ii, jj + d]
            A, b, r, off, m = A2, b2, r2, noff, mn
    json.dump(dict(modell=MODELL, N=N, P=P, werte=erg), open(HIER / f"v51_E4_{MODELL}_N{N}.json", "w", encoding="utf-8"))
    print(f"fertig ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()

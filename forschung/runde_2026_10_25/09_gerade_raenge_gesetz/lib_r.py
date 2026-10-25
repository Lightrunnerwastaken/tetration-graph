"""Rang 6 (Hexation, Basis b) mit dem gp-Motor: Kneser-Hexation hexx ueber der regulaeren Pentation pentx.
Fixpunkt L5 von pentx (Im > 0), Multiplikator mu = pent'(L5), Schroeder-Reihe sigma an L5,
regulaere Abel-Funktion h6 = log(sigma)/log(mu) und theta6(z) = h6(hex z) - z."""
import mpmath as mp
import motor

class RangR:
    def __init__(self, basis, dps, start, r, K=40, rho=mp.mpf("0.2"), rad=mp.mpf("0.04")):
        self.fv = "hyper(%d,%%s)" % (r - 1); self.fi = "ihyper(%d,%%s)" % (r - 1)
        self.B, self.dps = basis, dps
        mp.mp.dps = dps + 10
        self.K = K
        L = mp.mpc(start)
        for it in range(6):                      # Newton mit Taylor-Reihe um die Naeherung
            c = self.taylor(self.fv, L, rho, 128)
            # Nullstelle von f(L+u) - (L+u)
            u = mp.mpc(0)
            for _ in range(50):
                f = mp.polyval(c[::-1], u) - L - u
                d = mp.polyval([k * c[k] for k in range(len(c) - 1, 0, -1)], u) - 1
                du = f / d; u -= du
                if abs(du) < mp.mpf(10) ** (-dps - 5): break
            L = L + u
            if abs(u) < mp.mpf(10) ** (-dps + 3): break
        self.L = L
        c = self.taylor(self.fv, L, rho, 128)
        self.fk = [mp.mpc(0)] + c[1:K + 1]       # f(L+u) - L
        self.mu = self.fk[1]
        self.logmu = mp.log(self.mu)
        self.S = self._schroeder()
        self.rho = rho; self.rad = rad

    def taylor(self, fn, z0, r, M):
        pts = [z0 + r * mp.expjpi(2 * mp.mpf(j) / M) for j in range(M)]
        v = motor.werte(self.B, [fn % motor.gpzahl(p, self.dps + 10) for p in pts], self.dps)
        return [mp.fsum(v[j] * mp.expjpi(-2 * mp.mpf(j * k) / M) for j in range(M)) / M / r ** k for k in range(M // 2)]

    def _schroeder(self):
        K, f, mu = self.K, self.fk, self.mu
        def mul(a, b):
            c = [mp.mpc(0)] * (K + 1)
            for i in range(K + 1):
                if a[i] != 0:
                    for j in range(K + 1 - i): c[i + j] += a[i] * b[j]
            return c
        fp = [None, f[:K + 1]]
        for _ in range(2, K + 1): fp.append(mul(fp[-1], f[:K + 1]))
        s = [mp.mpc(0), mp.mpc(1)] + [mp.mpc(0)] * (K - 1)
        for m in range(2, K + 1):
            s[m] = mp.fsum(s[k] * fp[k][m] for k in range(1, m)) / (mu - mu ** m)
        return s

    def sigma_lokal(self, w):
        return mp.polyval(self.S[::-1], w - self.L)

    def h6(self, ws, rad=None, nmax=200):
        rad = self.rad if rad is None else rad
        """Abel-Funktion an L5 fuer eine Liste von Punkten (anziehend: vorwaerts mit pentx; abstossend: ipentx)."""
        vorw = abs(self.mu) < 1
        fn = self.fv if vorw else self.fi
        ws = list(ws); n = [0] * len(ws)
        for runde in range(nmax):
            offen = [i for i, w in enumerate(ws) if abs(w - self.L) > rad]
            if not offen: break
            neu = motor.werte(self.B, [fn % motor.gpzahl(ws[i], self.dps + 10) for i in offen], self.dps, speichern=False)
            for i, v in zip(offen, neu):
                ws[i] = v; n[i] += 1
        else:
            raise RuntimeError("h6: keine Konvergenz bei %d Punkten" % len(offen))
        s = -1 if vorw else 1
        return [mp.log(self.sigma_lokal(w)) / self.logmu + s * k for w, k in zip(ws, n)], n

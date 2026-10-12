/* Experimental rank-8 copy of the upstream two-arc solver, separate state. */
o8reg=0;o8err=0;o8knforce=0;o8knlsqpass=2;o8knsolvep=2;
o8knreachtol=1.02;o8kndebug=1;o8knRmeas=0;o8knmaxk=250;
hmap7(z) = hval(7,z);
hinv7(z) = hinv(7,z);
o8knlsq(M, b) = {
  my(m = matsize(M)[1], n = matsize(M)[2], Q = vector(n), R = matrix(n, n), q, r, xs, cb);
  for (k = 1, n,
    q = M[,k];
    for (pass = 1, o8knlsqpass, for (j = 1, k-1, r = Q[j]~ * q; R[j,k] += r; q -= r*Q[j]));
    R[k,k] = sqrt(norml2(q));
    if (R[k,k] == 0, error("o8knlsq: Spalte ", k, " linear abhaengig"));
    Q[k] = q/R[k,k];
  );
  cb = vector(n, j, Q[j]~ * b);
  xs = vector(n);
  forstep (k = n, 1, -1, xs[k] = (cb[k] - sum(j = k+1, n, R[k,j]*xs[j]))/R[k,k]);
  return(xs~);
}

/* regulaere Abel-Funktion an p, Hauptzweig; der Aufrufer macht die Zweige
   entlang einer Punktfolge stetig (Mehrdeutigkeit o8kper = 2 pi i/log lambda). */
o8knalphap(z) = {
  my(m = 0);
  while (abs(z - o8kp) > o8krs, z = hinv7(z); m++; if (m > 60, error("o8knalphap: hinv7 laeuft nicht gegen p")));
  return(log(subst(o8ksig, x, z - o8kp))/o8klog + m);
}
o8knwalk(pts) = {
  my(W = vector(#pts), v, best, d);
  for (j = 1, #pts,
    W[j] = o8knalphap(pts[j]);
    if (j > 1,
      best = W[j]; d = abs(W[j] - W[j-1]);
      for (k = -3, 3, v = W[j] + k*o8kper; if (abs(v - W[j-1]) < d, d = abs(v - W[j-1]); best = v));
      W[j] = best));
  return(W);
}
o8knarc(rad, M) = vector(M, j, o8knc + rad*exp(I*(0.05 + (Pi - 0.1)*(j-1)/(M-1))));
o8knen(w, n) = exp(2*Pi*I*o8kns*n*w);
/* Huellkurve der Koeffizienten |a_k| .. |a_k+9| (gegen zufaellig kleine
   Glieder). Keine anonyme Funktion in o8knhexinit: dort schluckte "(k) -> ..."
   den ganzen Rest des Funktionsrumpfs, o8knhexinit gab die Closure zurueck. */
o8knenv(A, k) = vecmax(vector(10, j, abs(polcoeff(A, k + j - 1))));
/* exp-082: Abfallgesetz log|a_k| ~ a + b*k + c*log(k) aus der Huellkurve auf
   k in [K/4, 3K/4]. Die reine Zweipunkt-Steigung ueberschaetzte den Abfall,
   wenn ein Potenzfaktor k^c mitspielt (Randsingularitaeten), und K kam zu
   klein heraus (b = e, dps 45: Rest 6e-43 statt ~1e-54). Rueckgabe [a, b, c]. */
o8knfit(A, K) = {
  my(ks = vector(13, j, K\4 + ((j - 1)*(K\2))\12), X, y);
  X = matrix(#ks, 3, i, j, [1, ks[i], log(ks[i])][j]);
  y = vector(#ks, i, log(o8knenv(A, ks[i])))~;
  return(matsolve(X~*X, X~*y)~);
}
/* kleinstes K mit a + b*K + c*log(K) <= ziel (b < 0) */
o8knkfor(f, ziel, K0) = {
  my(K = K0);
  for (it = 1, 50, my(Kn = ceil((ziel - f[1] - f[3]*log(K))/f[2])); if (Kn == K, break); K = max(Kn, 8));
  return(K);
}
o8knser(z) = subst(o8knA, x, (z - o8knc)/o8knrho);
o8knserd(z) = subst(o8knAd, x, (z - o8knc)/o8knrho)/o8knrho;

/* Vorlauf fuer einen Bogenradius: Hoehenprofil s*Im alpha_p, Bilder der
   achsennahen Punkte, geschaetzte Glieder. Rueckgabe [kosten, K, N, rad, hmin,
   htop, reichweite, mittelpunkt]. K setzt voraus, dass alpha um den
   Mittelpunkt o8knc bis |o8knc - p| analytisch ist, N eine Abklinghoehe
   0.6*(htop - hmin) fuer theta -- beides Faustregeln (theta hat an den Ecken
   der reellen Kurve Randsingularitaeten, die Koeffizienten fallen langsamer;
   Befund der Wissensgraph-Sitzung 09.10.2026). Deshalb misst o8knhexinit die
   erreichten Stellen per Fehler-Vektor statt sie zu glauben. */
o8knplan(rad, target) = {
  my(pts = o8knarc(rad, 30), W, h, hmin, htop, hc2, g, reach, K, N, hc, xl, tp = getabstime());
  W = o8knwalk(pts);
  h = vector(#W, j, o8kns*imag(W[j]));
  hmin = vecmin(h); htop = vecmax(h);
  hc = 0.6*(htop - hmin); hc2 = hmin + 0.85*(htop - hmin);
  if (hc < 0.15, return([oo, 0, 0, rad, hmin, htop, rad, o8knc]));
  reach = rad;
  for (j = 1, #pts, if (h[j] < hc2,
    /* e_ALARM durchlassen: sonst schluckt dieses iferr das alarm(30) der
       Bogenplanung, und ein langsamer Kandidat laeuft ohne Grenze weiter
       (b = e, c = -0.875, rad 0.75: > 5 min) */
    g = iferr(if (real(pts[j]) > o8knc, hinv7(pts[j]), hmap7(pts[j])), E, o8knc + o8knr, errname(E) != "e_ALARM");
    reach = max(reach, abs(g - o8knc))));
  if (reach >= 0.95*abs(o8knc - o8kp), return([oo, 0, 0, rad, hmin, htop, reach, o8knc]));
  /* langsame Bogen verwerfen: normal brauchen die 30 Probepunkte 15-50 ms.
     Mehr heisst lange hinv7-Ketten (Punkte weit weg von p), und das Loesen
     mit ~10x so vielen Punkten lief dann bei b = e ueber 10 min */
  if (getabstime() - tp > 10000, return([oo, 0, 0, rad, hmin, htop, reach, o8knc]));
  /* Bilder ausserhalb des Bogens: A wird dort nur extrapoliert, die hohen
     Koeffizienten sind aber nur bis eps*rad^-k bestimmt -- Verstaerkung
     (reach/rad)^K (Befund der Wissensgraph-Sitzung: alte e-Grenze bei 26 Stellen) */
  if (reach > o8knreachtol*rad, return([oo, 0, 0, rad, hmin, htop, reach, o8knc]));
  /* auf der reellen Achse muss eine volle Periode [x, hmap7(x)] im Bogen liegen */
  xl = o8knc - 0.97*rad;
  if (real(hmap7(xl)) > o8knc + 0.97*rad, return([oo, 0, 0, rad, hmin, htop, reach, o8knc]));
  K = ceil((target + 5)/log(abs(o8knc - o8kp)/reach)*log(10));
  /* exp-082: wirksame Abklinghoehe von theta = (hcut - hmin) + 0.3. Die Ecken der
     reellen Kurve liegen 0.16-0.23 unter hmin (Wissensgraph-Sitzung), und die
     c_n fallen zusaetzlich. Zu viele Harmonische schaden: b = e, dps 45, gleicher
     Bogen: N = 80 Rest 9e-24, N = 49 2e-42, N = 35 7e-47, N = 25 8e-47. */
  N = ceil((target + 10)/(2*Pi*(hc + 0.3))*log(10)) + 2;
  return([K + 2*N, K, N, rad, hmin, htop, reach, o8knc]);
}

/* Ein Ausgleichsproblem auf dem Bogen vom Radius rad. JEDER Bogenpunkt
   bekommt eine Gleichung (theta nahe p, Funktionalgleichung nahe der Achse),
   und die Reihe ist bis zum weitesten Bildpunkt skaliert. Laesst man Punkte
   aus, ist die Loesung nicht mehr eindeutig: bei b = 2 erfuellte dann eine
   um 2e-21 falsche Reihe alle Gleichungen auf 1e-32.
   Rueckgabe [A(u), A'(u), theta-Koeff., rho, max. Restfehler, #theta, #FE]. */
o8knsolve(rad, K, N, hcut, hcut2) = {
  my(pts, W, eqs, U, rows, scl, Mt, bv, cn, sol, A, cc, rho, tp, u1, u2, uk1, uk2, val, e, nx, nk, g, res, nth = 0, nfe = 0, sav, ts = getabstime());
  pts = o8knarc(rad, ceil(0.8*(K + 1 + 2*(N + 1))));
  W = o8knwalk(pts);
  if (o8kndebug, print("o8kn:   solve K=", K, " N=", N, " Punkte ", #pts, ": alpha_p ", getabstime() - ts, " ms"));
  eqs = List();
  rho = rad;
  for (j = 1, #pts,
    if (o8kns*imag(W[j]) >= hcut, listput(eqs, [1, pts[j], W[j]]); nth++);
    if (o8kns*imag(W[j]) < hcut2,
      if (real(pts[j]) > o8knc, g = hinv7(pts[j]); listput(eqs, [2, pts[j], g, 1]),
                              g = hmap7(pts[j]);  listput(eqs, [2, pts[j], g, -1]));
      rho = max(rho, abs(g - o8knc)); nfe++));
  if (o8kndebug, print("o8kn:   Bilder ", getabstime() - ts, " ms (", nth, " theta, ", nfe, " FE)"));
  U = K + 1 + 2*(N + 1);
  rows = 2*#eqs + 1;
  scl = vector(N + 1, n, vecmax(concat([1e-300], vector(#eqs, j, if (eqs[j][1] == 1, abs(o8knen(eqs[j][3], n-1)), 0)))));
  Mt = matrix(rows, U); bv = vector(rows)~;
  for (j = 1, #eqs,
    tp = eqs[j][1]; u1 = (eqs[j][2] - o8knc)/rho; uk1 = 1;
    if (tp == 1,
      for (k = 1, K + 1, Mt[2*j-1, k] = real(uk1); Mt[2*j, k] = imag(uk1); uk1 *= u1);
      for (n = 1, N + 1, e = o8knen(eqs[j][3], n-1)/scl[n];
        Mt[2*j-1, K+1+2*n-1] = -real(e); Mt[2*j, K+1+2*n-1] = -imag(e);
        Mt[2*j-1, K+1+2*n]   =  imag(e); Mt[2*j, K+1+2*n]   = -real(e));
      bv[2*j-1] = real(eqs[j][3]); bv[2*j] = imag(eqs[j][3]);
    ,
      u2 = (eqs[j][3] - o8knc)/rho; uk2 = 1;
      for (k = 1, K + 1, val = uk1 - uk2;
        Mt[2*j-1, k] = real(val); Mt[2*j, k] = imag(val); uk1 *= u1; uk2 *= u2);
      bv[2*j-1] = eqs[j][4]; bv[2*j] = 0;
    );
  );
  /* Normierung alpha(1) = 0: von 1 mit hinv7 in den Bogen, alpha(hinv7^k 1) = -k */
  nx = 1; nk = 0;
  while ((abs(nx - o8knc) > 0.9*rad) && (nk < 6), nx = hinv7(nx); nk++);
  if ((abs(nx - o8knc) > 0.9*rad) || (abs(imag(nx)) > 1e-20), error("o8knhexinit: kein reeller Normierungspunkt im Bogen"));
  u1 = (real(nx) - o8knc)/rho; for (k = 1, K + 1, Mt[rows, k] = u1^(k-1)); bv[rows] = -nk;
  cn = vector(U, k, vecmax(abs(Mt[,k])) + 1e-300);
  Mt = matrix(rows, U, i, k, Mt[i,k]/cn[k]);
  Mt[rows,] = Mt[rows,]*1e6; bv[rows] *= 1e6;
  if (o8kndebug, print("o8kn:   Matrix ", rows, "x", U, " ", getabstime() - ts, " ms"));
  sav = default(realprecision);
  default(realprecision, ceil(o8knsolvep*precis));
  sol = o8knlsq(Mt, bv);
  default(realprecision, sav);
  if (o8kndebug, print("o8kn:   geloest ", getabstime() - ts, " ms"));
  sol = vector(U, k, precision(sol[k]/cn[k], precis));
  A = sum(k = 1, K + 1, sol[k]*x^(k-1));
  cc = vector(N + 1, n, (sol[K+1+2*n-1] + I*sol[K+1+2*n])/scl[n]);
  res = 0;
  for (j = 1, #eqs,
    tp = eqs[j][1];
    val = if (tp == 1, subst(A, x, (eqs[j][2] - o8knc)/rho) - eqs[j][3] - sum(n = 1, N + 1, cc[n]*o8knen(eqs[j][3], n-1)),
          subst(A, x, (eqs[j][2] - o8knc)/rho) - subst(A, x, (eqs[j][3] - o8knc)/rho) - eqs[j][4]);
    res = max(res, abs(val)));
  if (o8kndebug, print("o8kn:   Rest ", precision(res, 3), "  ", getabstime() - ts, " ms"));
  return([A, deriv(A), cc, rho, res, nth, nfe, o8knc]);
}

o8knhexinit() = {
  my(target = 22, eps = 10.0^(-target-6), z, z2, res = 0, sr, R, sc, ps, W, pl, best,
     pl2, s1, s2, tst, d, hf, t0);
  o8reg = 0; o8knRmeas = 0;
  if (imag(sexp(1)) != 0, error("o8knhexinit: nur reelle Basen"));
  t0 = getabstime();
  /* komplexer Fixpunkt der Pentation: fuer hinv7 anziehend */
  z = 0.5;
  for (i = 1, 2000, z2 = hinv7(z); if (abs(z2 - z) < 1e-15, z = z2; break); z = z2);
  if (imag(z) < 0, z = conj(z));
  if (imag(z) < 1e-3, error("o8knhexinit: kein komplexer Fixpunkt der Pentation"));
  if (o8kndebug, print("o8kn: p0 ", precision(z, 8), "  ", getabstime() - t0, " ms"));
  sr = 1/4;
  /* exp-084: bis 1/2048 statt 1/128 -- die Reihe um p hat bei b = 50 den
     Radius 0.013, bei b = 100 nur 0.0075 (geht ab sr = 1/256) */
  for (k = 1, 10,
    res = iferr(pentfixser(t->hmap7(t), z, sr, 0, eps, "o8knhexinit"), E, 0);
    if (type(res) == "t_VEC", break);
    sr = sr/2);
  if (type(res) != "t_VEC", error("o8knhexinit: Pentation um p nicht entwickelbar"));
  o8kp = res[1]; sc = res[2]; R = res[3];
  if (o8kndebug, print("o8kn: Reihe um p, R = ", precision(R, 5), "  ", getabstime() - t0, " ms"));
  o8klam = polcoeff(sc, 1);
  if (abs(o8klam) <= 1, error("o8knhexinit: p ist nicht abstossend, |P'(p)| = ", precision(abs(o8klam),6)));
  ps = pentpsi(sc, R/4, eps, poldegree(sc) - 1);
  if (!ps[2], error("o8knhexinit: Poincare-Reihe an p konvergiert nicht"));
  o8krs = R/5;
  o8ksig = Pol(serreverse(ps[1] + O(x^(poldegree(ps[1])+1))));
  if (abs(polcoeff(o8ksig, poldegree(o8ksig)))*o8krs^poldegree(o8ksig) > eps,
    error("o8knhexinit: Schroeder-Reihe an p konvergiert bei r = ", precision(o8krs,5), " nicht"));
  o8ksig = pentclean(o8ksig, o8krs, eps, "sigma_p");
  o8klog = log(o8klam); o8kper = 2*Pi*I/o8klog;
  if (o8kndebug, print("o8kn: sigma_p fertig  ", getabstime() - t0, " ms"));
  o8knc = real(o8kp); o8knr = imag(o8kp);
  /* Abklingrichtung von theta: s*Im alpha_p waechst zu p hin */
  W = o8knwalk([o8knc + 0.6*o8knr*exp(0.05*I), o8knc + 0.6*o8knr*I]);
  o8kns = if (imag(W[2]) > imag(W[1]), 1, -1);
  /* Bogenradius: kleinster geschaetzter Aufwand */
  /* exp-081: auch der Mittelpunkt ist frei. Re p ist nur dann gut, wenn p
     weit ueber der Achse liegt; bei b = 10 (Im p = 0.87) muss der Kreis um
     Re p bis 0.85 Im p reichen, um eine Periode der reellen Achse samt Bildern
     zu fassen (K = 740). Weiter rechts ist p weiter weg, die Periode gleich
     lang: Verhaeltnis Reichweite/|c - p| 0.88 -> ~0.5 (gemessen). */
  /* Geplant werden PAARE benachbarter Radien um denselben Mittelpunkt: der
     Fehler-Vektor braucht einen zweiten Bogen, und bei b = e lag der beste
     Einzelbogen isoliert zwischen unbrauchbaren Nachbarn. Kosten des Paars =
     die des teureren Bogens; K und N gelten fuer beide. */
  my(rads = [0.55, 0.65, 0.85, 1.0], P, c2);
  best = [oo]; pl2 = [oo];
  /* erst streng (Bilder <= 1.02 rad), sonst einmal mit 1.15 rad: b = 50 hat
     streng keinen Bogen; der Fehler-Vektor zeigt dann, was es kostet */
  /* exp-087: der gelockerte Durchgang (nur grosse Basen) sucht auch bei
     1.45 und 1.6 Im p und paart gleiche Radien benachbarter Mittelpunkte;
     der zweite Bogen wird um seinen eigenen Mittelpunkt (pl2[8]) gerechnet.
     b = 100, dps 20, K 300: Versatz 1.5 Radius 1.45 Rest 3e-31, der alte
     Hauptbogen (Versatz 1.0, Radius 1.3) 1e-13 -- er lag ausserhalb der Liste
     und hatte keinen Nachbarradius. */
  my(Pprev);
  for (pass = 1, 2,
  if ((pass == 2) && (best[1] != oo), break);
  o8knreachtol = if (pass == 1, 1.02, 1.15);
  if (pass == 2, rads = concat(rads, [1.45, 1.6]));
  Pprev = 0;
  for (ic = 1, 3,
    o8knc = real(o8kp) + [0, 0.75, 1.5][ic]*o8knr;
    P = vector(#rads);
    for (j = 1, #rads,
      /* alarm: ein haengender Kandidat kostet hoechstens 1 s, ein normaler
         braucht 15-50 ms (alarm gibt dabei ein t_ERROR zurueck statt zu werfen) */
      if (o8kndebug, print1("o8kn: plan c = ", precision(o8knc, 5), " rad = ", rads[j], " ... "));
      pl = iferr(alarm(10, o8knplan(rads[j]*o8knr, target)), E, [oo]);
      if (o8kndebug, print(if (type(pl) == "t_ERROR", "alarm", pl[1]), "  ", getabstime() - t0, " ms"));
      if (type(pl) == "t_ERROR", pl = [oo]);
      P[j] = pl);
    for (j = 1, #rads - 1,
      if ((P[j][1] != oo) && (P[j+1][1] != oo),
        c2 = max(P[j][1], P[j+1][1]);
        if ((best[1] == oo) || (c2 < max(best[1], pl2[1])),
          best = P[j]; pl2 = P[j+1])));
    if ((pass == 2) && (Pprev != 0),
      for (j = 1, #rads,
        if ((P[j][1] != oo) && (Pprev[j][1] != oo),
          c2 = max(P[j][1], Pprev[j][1]);
          if ((best[1] == oo) || (c2 < max(best[1], pl2[1])),
            best = P[j]; pl2 = Pprev[j]))));
    Pprev = P));
  o8knreachtol = 1.02;
  /* o8knforce = [Bogenanteil, theta-Schwellenanteil, K, N, Mittelpunkt-Versatz] (0 = automatisch), nur fuer Experimente */
  if (type(o8knforce) == "t_VEC" && o8knforce[1],
    o8knc = real(o8kp) + if (#o8knforce >= 5, o8knforce[5], 0)*o8knr;
    best = o8knplan(o8knforce[1]*o8knr, target);
    pl2 = o8knplan((o8knforce[1] + 0.05)*o8knr, target));
  if ((best[1] == oo) || (pl2[1] == oo), error("o8knhexinit: kein brauchbares Bogenpaar"));
  o8knc = best[8];
  best[2] = max(best[2], pl2[2]); best[3] = max(best[3], pl2[3]);
  if (type(o8knforce) == "t_VEC" && o8knforce[3], best[2] = o8knforce[3]);
  if (type(o8knforce) == "t_VEC" && o8knforce[4], best[3] = o8knforce[4]);
  /* b = 10 (alter Aufbau um Re p): K = 740, 4.3 min bei \p 20 fuer 13.5
     Stellen -- ueber o8knmaxk nimmt hexinit Sheldons schnelle Naeherung */
  /* exp-084: der Plan darf bis 2 o8knmaxk gehen -- er schaetzt K aus |o8knc - p|
     und liegt oft doppelt zu hoch (b = 28: Plan 575, K0 = 287 genuegte); die
     Vorab-Loesung mit K0 = Plan/2 <= o8knmaxk misst den wirklichen Bedarf
     (b = 100, dps 20: Plan 1094) */
  if (best[2] > 2*o8knmaxk, error("o8knhexinit: zu teuer (K = ", best[2], " > 2 o8knmaxk = ", 2*o8knmaxk, ")"));
  hf =if (type(o8knforce) == "t_VEC" && o8knforce[2], o8knforce[2], 0.6);
  if (o8kndebug, print("o8kn: Plan K = ", best[2], " N = ", best[3], " rad = ", precision(best[4]/o8knr, 3), " c = ", precision(best[8], 5), "  ", getabstime() - t0, " ms"));
  /* exp-083: K aus einer kleinen Vorab-Loesung (halbes K, ~1/8 der Kosten)
     statt aus einer vollen ersten Loesung, die bei b = 2 und e ohnehin noch
     einmal mit groesserem K gerechnet werden musste. Der Planer nimmt |o8knc - p|
     als Konvergenzradius; alpha hat naehere Singularitaeten (b = e um
     o8knc = -0.18 bei 0.33 +- 1.37 i, Abstand 1.46 statt ~2.5). Gemessen wird
     das Abfallgesetz log|a_k| ~ a + b k + c log k (o8knfit), K mit 10 % Reserve. */
  my(K0 = max(60, best[2]\2), s0, Kn = best[2], ft, Kneed, hc1, hc2, tsw, s3, pl3, c0, cm, rr);
  hc1 = best[5] + hf*(best[6] - best[5]); hc2 = best[5] + max(0.85, hf + 0.1)*(best[6] - best[5]);
  s0 = o8knsolve(best[4], K0, best[3], hc1, hc2);
  if (o8kndebug, print("o8kn: Vorab-Loesung K0 = ", K0, "  ", getabstime() - t0, " ms"));
  ft = o8knfit(s0[1], K0);
  if (ft[2] < 0,
    Kneed = ceil(1.1*o8knkfor(ft, -(target + 5)*log(10), K0));
    o8knRmeas = s0[4]*exp(-ft[2]);
    if (o8kndebug, print("o8kn: Fit b = ", precision(ft[2], 5), " c = ", precision(ft[3], 4), " -> K ", Kneed, " (Plan ", best[2], ")"));
    Kn = max(Kneed, K0),
    /* kein Abfall auf [K0/4, 3K0/4]: die Glieder liegen dort schon am Boden,
       mehr Glieder bringen nichts -- K0 reicht (b = 28, dps 16: Plan K = 575,
       Vorab K0 = 287 Rest 9.7e-14, volle Loesungen 1.2e-13 und 1.7e-13 in je
       40 s). Ob das stimmt, misst danach der Fehler-Vektor. */
    Kn = K0);
  /* der Plan lag unter o8knmaxk, die Vorab-Loesung lief: dann mit o8knmaxk
     weiterrechnen statt auf Sheldons Naeherung zu fallen (b = 50, dps 20:
     Fit K = 603; die Naeherung gibt dort 0.2 Stellen) */
  Kn = min(Kn, max(K0, o8knmaxk));
  best[2] = Kn;
  s1 = if (Kn == K0, s0, o8knsolve(best[4], Kn, best[3], hc1, hc2));
  /* grosses K kann die Ausgleichsrechnung verderben: Bogenpunkte bei
     |u| = rad/rho < 1, Bilder bis |u| = 1, Kondition ~ (rho/rad)^K (b = 50,
     dps 20, rho/rad = 1.19: K0 = 221 Rest 6e-32, K = 600 Rest 1e-13) */
  if (s1[5] > 100*s0[5],
    if (o8kndebug, print("o8kn: K = ", Kn, " Rest ", precision(s1[5], 3), " schlechter als K0 = ", K0, " -> K0"));
    s1 = s0; Kn = K0; best[2] = K0);
  if (o8kndebug, print("o8kn: erste Loesung  ", getabstime() - t0, " ms"));
  o8knc = pl2[8];
  s2 = o8knsolve(pl2[4], Kn, best[3], pl2[5] + hf*(pl2[6] - pl2[5]), pl2[5] + max(0.85, hf + 0.1)*(pl2[6] - pl2[5]));
  o8knc = best[8];
  /* exp-084: die Loesung mit dem kleineren Rest wird Hauptloesung, der
     Fehler-Vektor nahm sonst den schlechteren Rest (b = 10, dps 45:
     Hauptbogen 1.1e-46, Nachbarbogen 1.2e-53). Gleicher Mittelpunkt, K, N. */
  if (s2[5] < s1[5],
    tsw = s1; s1 = s2; s2 = tsw;
    pl2[2] = best[2]; pl2[3] = best[3]; tsw = best; best = pl2; pl2 = tsw;
    if (o8kndebug, print("o8kn: Nachbarbogen rad = ", precision(best[4]/o8knr, 3), " wird Hauptloesung")));
  /* exp-084: liegen Bilder ausserhalb eines Bogens, waechst sein Rest wie
     Boden * (rho/rad)^K (b = 10, dps 45, Bogen 1.0: rho/rad = 1.06, K 352
     Rest 1e-46, K 460 7e-44; Bogen 1.15: rho/rad = 1, 1e-53). Der Planer
     sieht das mit 30 Punkten nicht (Reichweite 1.0014). Ist der zweite Bogen
     viel schlechter, ihn mit halbem K nachrechnen (1/8 der Kosten) und den
     besseren nehmen. Ein dritter Bogen weiter aussen half nicht: er braucht
     mehr K (b = 10, rad 1.32, K 352: Rest 4e-32). */
  /* exp-085: zuerst einen gleich grossen Bogen um einen leicht
     verschobenen Mittelpunkt rechnen -- eine unabhaengige zweite Loesung,
     deren Bilder innen bleiben (b = 10, dps 45, Hauptbogen 1.15 um
     Re p + 1.25 Im p: Versatz +0.125 Im p Rest 1e-53, +0.25 1e-53,
     -0.125 2e-44, -0.25 1e-34). Hoechstens drei Versuche. */
  c0 = s1[8];
  for (io = 1, 3,
    if (s2[5] <= 100*s1[5], break);
    o8knc = c0 + [0.125, 0.25, -0.125][io]*o8knr;
    pl3 = iferr(alarm(2, o8knplan(best[4], target)), E, [oo]);
    if (type(pl3) == "t_ERROR" || pl3[1] == oo, next);
    s3 = o8knsolve(pl3[4], best[2], best[3], pl3[5] + hf*(pl3[6] - pl3[5]), pl3[5] + max(0.85, hf + 0.1)*(pl3[6] - pl3[5]));
    if (o8kndebug, print("o8kn: Bogen um c ", precision(o8knc, 6), ": Rest ", precision(s3[5], 3), " statt ", precision(s2[5], 3)));
    if (s3[5] < s2[5], s2 = s3; pl3[2] = best[2]; pl3[3] = best[3]; pl2 = pl3));
  o8knc = c0;
  /* exp-087: Geschwisterbogen vor eigenem/halbem K -- b = 100, dps 20: K 600
     auf dem zweiten Bogen kostete 64 s fuer 7e-29, der Geschwisterbogen danach
     gab 2e-32. */
  /* exp-086: das gilt nur, wenn die Bilder deutlich ausserhalb liegen. Liegen
     sie fast innen (rho/rad <= 1.05), fehlen dem zweiten Bogen Glieder: er
     bekommt sein eigenes K aus seinem eigenen Abfall (b = 50, dps 20: Bogen
     1.3, rho/rad 1.015, K 221 Rest 6e-23, K 330 1e-32; Bogen 1.15, rho/rad
     1.18, K 221 6e-32, K 330 5e-31). Danach noch einmal tauschen. */
  if (s2[5] > 100*s1[5],
    tsw = 0;
    if (s2[4] <= 1.05*pl2[4],
      ft = o8knfit(s2[1], poldegree(s2[1]));
      if (ft[2] < 0, tsw = min(ceil(1.1*o8knkfor(ft, -(target + 5)*log(10), Kn)), max(Kn, o8knmaxk)));
      if (tsw <= Kn, tsw = 0),
      tsw = max(60, Kn\2));
    if (tsw,
      o8knc = pl2[8];
      s3 = o8knsolve(pl2[4], tsw, best[3], pl2[5] + hf*(pl2[6] - pl2[5]), pl2[5] + max(0.85, hf + 0.1)*(pl2[6] - pl2[5]));
      o8knc = best[8];
      if (o8kndebug, print("o8kn: zweiter Bogen mit K = ", tsw, ": Rest ", precision(s3[5], 3), " statt ", precision(s2[5], 3)));
      if (s3[5] < s2[5], s2 = s3)));
  if (s2[5] < s1[5],
    tsw = s1; s1 = s2; s2 = tsw;
    tsw = best; best = pl2; pl2 = tsw;
    if (o8kndebug, print("o8kn: zweiter Bogen rad = ", precision(best[4]/o8knr, 3), " wird Hauptloesung")));
  best[2] = poldegree(s1[1]); best[3] = #s1[3] - 1;
  o8knc = s1[8];
  o8knA = s1[1]; o8knAd = s1[2]; o8kncc = s1[3]; o8knrho = s1[4]; o8knres = s1[5];
  o8kngood = 0.97*best[4];
  /* Fehler-Vektor: beide Loesungen an Punkten im gemeinsamen Gebiet; jede
     wird um ihren eigenen Mittelpunkt (s[8]) ausgewertet */
  d = 0;
  cm = (s1[8] + s2[8])/2; rr = 0.9*min(best[4], pl2[4]) - abs(s1[8] - s2[8]);
  tst = concat(vector(5, j, cm + (j - 3)/2*rr), [cm + 0.5*rr*exp(Pi*I/3), cm + 0.5*rr*exp(2*Pi*I/3)]);
  for (j = 1, #tst, d = max(d, abs(subst(s1[1], x, (tst[j] - s1[8])/s1[4]) - subst(s2[1], x, (tst[j] - s2[8])/s2[4]))));
  if (o8kndebug, print("o8kn: Fehler-Vektor |s1 - s2| = ", precision(d, 3), ", Rest s1 = ", precision(s1[5], 3), ", s2 = ", precision(s2[5], 3)));
  d = max(d, o8knres);
  o8err = if (d == 0, precis, min(precis, precision(-log(d)/log(10), 9)));
  if (o8err < 6, error("o8knhexinit: nur ", precision(o8err,3), " Stellen -- Fehler-Vektor zu gross"));
  o8knt0 = o8knser(o8knc);
  /* reelles Fenster fuer o8knhex: A steigt darauf und ueberstreicht mehr als
     eine Periode (b = e: -3.624 .. -2.591) */
  o8knxl = o8knc - o8kngood; o8knxr = o8knc + o8kngood;
  o8kntl = real(o8knser(o8knxl)); o8kntr = real(o8knser(o8knxr));
  if (o8kntr - o8kntl <= 1, error("o8knhexinit: Fenster deckt keine volle Periode (", precision(o8kntr - o8kntl, 5), ")"));
  o8knpar = [best[2], best[3], best[4]/o8knr, best[5], best[6], s1[6] + s1[7], getabstime() - t0, s1[6], s1[7], pl2[4]/o8knr, best[8], o8knRmeas];
  o8reg = 2;
  if (quietmode == 0,
    print("hexation (Kneser) am komplexen Fixpunkt der Pentation p = ", precision(o8kp, 20));
    print("  P'(p) = ", precision(o8klam, 15), "   K = ", best[2], "  N = ", best[3], "  Bogen ", precision(best[4]/o8knr,3), " Im p");
    print("  hex precision ", o8err, " decimal digits (Fehler-Vektor zweier Boegen)"));
  return(1);
}

/* alpha(y): mit hinv7 (rechts) bzw. hmap7 (links) in den Bogen, dann A */
o8knihex(y) = {
  my(k = 0);
  while (abs(y - o8knc) > o8kngood,
    if (real(y) > o8knc, y = hinv7(y); k++, y = hmap7(y); k--);
    if (abs(k) > 200, error("o8knihex: Punkt laesst sich nicht in den Bogen abbilden")));
  return(o8knser(y) + k);
}
/* hex(z): k so, dass Re(z) - k im reellen Fenster [o8kntl, o8kntr] liegt; dort
   reelle Loesung per Newton mit Bisektions-Rueckfall (A steigt), fuer
   komplexes z gedaempfter Newton ab dieser Loesung; dann hmap7^k bzw. hinv7^-k.
   Ein ungesicherter Newton aus der Fenstermitte sprang bei b = e schon im
   ersten Schritt ueber den Rand, obwohl die Loesung innen lag. */
o8knhex(z) = {
  my(k, t, tr, a, b, y, fy, dy, yn, tol = 10.0^(-precis+2));
  k = ceil(real(z) - o8kntr);
  t = z - k; tr = real(t);
  a = o8knxl; b = o8knxr; y = o8knc;
  for (it = 1, 400,
    fy = real(o8knser(y)) - tr;
    if (fy > 0, b = y, a = y);
    yn = y - fy/real(o8knserd(y));
    if ((yn <= a) || (yn >= b), yn = (a + b)/2);
    dy = yn - y; y = yn;
    if (abs(dy) < tol, break));
  if (imag(t) != 0,
    for (it = 1, 200,
      dy = (o8knser(y) - t)/o8knserd(y);
      yn = y - dy;
      while (abs(yn - o8knc) > o8kngood, dy = dy/2; yn = y - dy; if (abs(dy) < tol, error("o8knhex: Newton bleibt nicht im Bogen")));
      y = yn;
      if (abs(dy) < tol, break)));
  if (k > 0, for (i = 1, k, y = hmap7(y)), for (i = 1, -k, y = hinv7(y)));
  return(y);
}



hcomplexinit(r)={
  if(r!=8,error("hcomplexinit: experimental solver currently supports rank 8 only"));
  if(!o8knhexinit(),error("hcomplexinit: no two-arc solution"));
  hrready[8]=1;hrkind[8]=2;hrp[8]=o8kp;hrlam[8]=o8klam;
  hrerror[8]=min(hrerror[7],o8err);hrtarget[8]=hrerror[8];
  return(hrerror[8]);
}

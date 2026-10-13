# Recherche: Wikis und Tetration-Forum gegen den Graphen (08.10.2026, Vollabgleich 10.10.2026)

## Durchsucht

| Quelle | Umfang | Bemerkung |
|---|---|---|
| Wikipedia en: Tetration, Iterated logarithm (= Super-logarithm), Iterated function | vollständig | |
| Wikipedia de: Potenzturm, Hyper-Operator; Wikipedia ru: Тетрация | vollständig | |
| Googology Wiki: Tetration | vollständig | über den eingebauten Browser (Abruf sonst gesperrt) |
| MathWorld: Power Tower | vollständig | |
| Mizugadro/TORI (Kouznetsov): Tetration | Archivfassung Dez. 2024 | Seite selbst seit April 2026 offline |
| Citizendium: Tetration | **nicht gelesen** | Bot-Prüfung (Cloudflare) — nicht umgangen |
| Tetration-Forum: FAQ (Trappmann/Robbins 2008), Forums-Wiki (17 Seiten), Liste offener Probleme (TPID) | vollständig | |
| Tetration-Forum: Titel aller 1260 Themen (Diskussion, Hyperoperationen, Rechnen) | alle Titel | rund 60 Themen gelesen |

## Kurzfazit

1. **Vorarbeiten zu unseren eigenen Ergebnissen** — das Wichtigste. Mehrere Ergebnisse, die im Graphen als „vermutlich neu“ stehen, hat das Forum
   schon früher gesehen. Das muss in die Bewertungen, bevor die Bibliothek öffentlich wird.
2. **Etwa 25 belegte, prüfbare Aussagen**, die fehlen. Die meisten sind elementar oder klassisch: Superwurzeln, Pentation, Tetra-Euler-Zahl,
   Geschichte, Grenzfälle, Gelfond–Schneider-Folgerungen.
3. **Etwa 12 Vermutungen und Fragen** aus der TPID-Liste und den Foren-Diskussionen, die als Vermutung oder Frage gut hineinpassen.
4. **Ein neuer eigener Befund während der Recherche:** TPID 6 mit genauer Rate, $x_n-\eta\approx 2\pi^2\eta/(e\,n^2)$, erklärt durch unser Engpass-Gesetz.

---

## A. Vorarbeiten und Korrekturen an bestehenden Knoten (zuerst erledigen)

| Knoten | Was im Forum/Wiki steht | Quelle | Folge |
|---|---|---|---|
| `kl.beob.zwei_abel` | Dass die regulären Iterationen am unteren und oberen Fixpunkt von $\sqrt2^x$ sich um etwa $10^{-24}$ unterscheiden, fanden Asimov, Hickerson und Schroeppel um 1991 (unveröffentlicht); im Forum 2007 als „Bummer!“ nachgerechnet. | Forum tid 69 | Quelle ergänzen, Neuheit „bekannt“. Unser Beitrag: 170 Stellen, $\lvert\hat p(1)\rvert=3.63\cdot10^{-25}$, das Amplitudengesetz. |
| `bw.satz.phasengesetz`, `bw.prop.hub` | Jay Fox (2007): Basiswechsel-Formel ${}^xa=\lim\log_a^{\circ n}({}^{n+x+\mu_b(a)}b)$ mit **additiver** Konstante, $\mu_a(c)=\mu_a(b)+\mu_b(c)$. Sheldon Levenstein (2009): $\mathrm{slog}_2-\mathrm{slog}_e\approx1.1282$, $\mathrm{slog}_3-\mathrm{slog}_e\approx-0.1926$, $\mathrm{slog}_{10}-\mathrm{slog}_e\approx-1.1364$, und: der Wechsel „wackelt ein wenig“, ist also nicht konstant. | Forum tid 14, 236 | Quellen ergänzen. Mittel und Welle sind im Forum **beobachtet**; neu bei uns sind der Beweisrahmen (Schwänze, Kreisabbildung, Kozykel, Antwortkern) und die Präzision. Die Werte stimmen mit $m(2)$, $m(3)$, $m(10)$ überein. |
| `bw.verm.streifen` ↔ `bw.frage.streifen_gevrey` | Sheldon (2011, MSE und Forum) **vermutet**: Jede Tetration, die man per Basiswechsel aus einer analytischen Lösung einer anderen Basis gewinnt, ist $C^\infty$, aber **nirgends analytisch** — die Singularitäten rücken mit jeder Logarithmenstufe an die reelle Achse (Radius 0.46 → 0.035 → …). Trappmann (2009): $\log_a^{\circ n}\circ\exp_b^{\circ n}$ hat Verzweigungspunkte bei $\log_b^{\circ n}\{a,\dots,\exp_a^{\circ n-2}(1)\}$. Ungelöst (Stand Forum). | Forum tid 333, 578, 825, 1188 | Neue Gegen-Vermutung aufnehmen, die `bw.verm.streifen` widerspricht. Passt zu unserem Gevrey-Befund ($\sigma^*\to0$). Die Streifen-Vermutung als „umstritten“ markieren. |
| `kn.konstr.walker` | Walkers $C^\infty$-Lösung ist ein **Basiswechsel-Grenzwert** von $e^x-1$ zu $\exp$ ($\mathrm{dexp}^{-n}\circ\exp^n$); bewiesen ist nur $C^\infty$ auf der reellen Achse, vermutet: nirgends analytisch. | Forum tid 333, 825 | Text präzisieren, Verweis auf Basiswechsel. |
| `kl.satz.eta_koordinate` | Trappmann–Kouznetsov (Math. Comp. 2012) geben die Entwicklung $F(z)=e\bigl(1-\frac2z(1+\sum P_m(-\ln z)/(3z)^m)\bigr)$ der beiden regulären Superexponentiale zur Basis $\eta$ — dieselbe Struktur wie unsere Kehrwert-Koordinate. | Mizugadro | Quelle ergänzen, Neuheit „bekannt“ für die Entwicklung; neu bleiben bei uns die exakten rationalen $c_k$ und $C_T$ auf 30 Stellen. |
| `kl.satz.parabolisch` | Das **obere** Superexponential zur Basis $\eta$ („cheta“, Jay Fox 2009) ist ganz; $\check\eta(z)=\exp_\eta^{\circ z}(e^2)$ mit $\check\eta(-1)=2e$, $\check\eta(0)=e^2$, $\check\eta(1)=e^e$. | Forum tid 260, 316 | Eigene Aussage „obere Superexponentialfunktion“ (siehe B). |
| `kl.beob.stetig_eta` | Kouznetsov und mike3 halten $b=\eta$ für einen **Verzweigungspunkt** der Tetration als Funktion der Basis. Kouznetsovs vorläufige Taylor-Tabelle in $b-\eta$ hat $c_{0,1}=0.61061$ — passt zu unserem $T_\eta'(0)=0.611095$. | Mizugadro, Forum tid 648 | Als offene Frage mit Literatur verbinden: glatt ($C^\infty$), aber nicht analytisch in $b$? |
| `kl.satz.eta_digamma` | Gottfried Helms fragte 2010, ob $\sum(e-\eta\uparrow\uparrow k)$ konvergiert. | Forum tid 476 | Unser Satz beantwortet das: divergiert wie $2e\ln n$. Als „beantwortet Forumsfrage“ vermerken. |
| `kn.satz.hooshmand` | Wortlaut jetzt aus Wikipedia: $f$ stetig, $f(x)=a^{f(x-1)}$, $f(0)=1$, differenzierbar auf $(-1,0)$, $f'$ monoton dort, $f'(0^+)=\ln a\,f'(0^-)$ ⇒ linear. Zweiter Satz (Basis $e$): konvex auf $(-1,0)$ und $f'(0^-)\le f'(0^+)$ ⇒ linear. Folge: Jede andere hinreichend glatte Lösung hat einen **Wendepunkt in $(-1,0)$**. | Wikipedia en | „Wortlaut ungeprüft“ durch Sekundärzitat ersetzen; Wendepunkt-Folge mit `gr.beob.wendepunkte` verbinden. |
| `fp.satz.shell_thron`, `fp.satz.turm_reell` | Geschichte: Euler 1783, Eisenstein 1844, Barrow 1936 ($\lvert\log\lambda\rvert\le1$), Thron 1957, Shell 1962, Macintyre 1966 ($i^{i^{\cdots}}$), Knoebel 1981 (Übersicht). | MathWorld, Wikipedia de | Quellen ergänzen. |
| `kn.verm.carleman_kneser` | TPID 18 (Helms 2022): Konvergiert die „Polynom-Tetration“ aus abgeschnittenen Carleman-Matrizen gegen Kneser? — die sexp-Seite unserer slog-Vermutung. | Forum tid 162 | Verwandte Vermutung verlinken. |

## B. Neue Aussagen, die in die Bibliothek gehören (belegt und prüfbar)

Priorität: **1** = gehört zum Kern einer Bibliothek, **2** = gute Ergänzung, **3** = Randthema.

| Prio | Vorschlag (Gebiet) | Inhalt | Prüfung |
|---|---|---|---|
| 1 | Superwurzel (gr, neuer Strang „Umkehrungen“) | $\mathrm{ssrt}(x)=e^{W(\ln x)}=\ln x/W(\ln x)$; Anzahl reeller Lösungen (zwei für $e^{-1/e}<x<1$, eine für $x>1$, keine darunter); Reihe um $1$: $1+(x-1)-(x-1)^2+\tfrac32(x-1)^3-\dots$; dritte Superwurzel über $x_{n+1}=\exp W(W(x_n\ln y))$ | direkt nachrechnen |
| 1 | Superwurzel arithmetisch (ar) | Gelfond–Schneider: $\sqrt[s]{n}$ (zweite Superwurzel einer natürlichen Zahl) ist ganz oder transzendent; die dritte ist ganz oder irrational; offen, ob irrationale Superwurzeln transzendent sind | Beweis zwei Zeilen |
| 1 | Unendliche Superwurzel (fp) | ${}^\infty\!\sqrt{x}=x^{1/x}$ für $1/e\le x\le e$; der unendliche Turm ist streng wachsend und bijektiv $[e^{-e},\eta]\to[1/e,e]$ | elementar |
| 1 | **TPID 6 mit Rate** (fp/bw, neu) | $x_n$ mit ${}^nx_n=n$ geht gegen $\eta$ (im Forum bewiesen); **neu**: $x_n-\eta\approx2\pi^2\eta/(e\,n^2)=10.49/n^2$ — gemessen $1.05\cdot10^{-5}$ bei $n=10^3$, $1.05\cdot10^{-7}$ bei $n=10^4$; erklärt durch das Engpass-Gesetz | Experiment schreiben |
| 1 | Pentation (neuer Strang „Höhere Stufen“) | Reeller Fixpunkt von $T_e$ bei $-1.85035452902718\ldots$; Pentation zur Basis $e$ als reguläre Iteration dort (Fox, Robbins, Kouznetsov); $e\uparrow\uparrow\uparrow(-1)=0$, $(-2)\mapsto-1$; für $b<\eta$ zweiter Fixpunkt $>0$ | mit unserem Kneser-Motor |
| 1 | Tetra-Euler-Zahl und tetra-kritische Basis (Höhere Stufen) | $\eta_4=1.63532449671527\ldots$, $e_4=3.08853227180671\ldots$ mit drei Kennzeichnungen: Grenzwert unendlicher Pentation, Maximum der Selbst-Tetrawurzel, Übergang von drei zu einem Fixpunkt von $b\uparrow\uparrow x$ — wie $\eta$ und $e$ eine Stufe tiefer | mit Kneser nachrechnen |
| 1 | Obere Superexponentialfunktion (kl) | Für $1<b<\eta$ am oberen Fixpunkt, normiert $\mathrm{usexp}(0)=a+1$; bei $\eta$: cheta, ganz | mit `AmFixpunkt` |
| 1 | q-Binomial-Formel (kl) | Reguläre Tetration aus den ganzzahligen Werten: $T(z)=\sum_n\sum_k(-1)^{n+k}q^{\binom{n-k}2}\binom zn_q\binom nk_qT(k)$, $q=\lambda$ (Reshetnikov, MO 259278). **Geprüft**: $\sqrt2$, $z=\tfrac12$: Abstand $2\cdot10^{-16}$ bei $N=40$; danach wachsen die Glieder (Summation nötig) — *Korrektur Runde 2: Rechenartefakt; mit genug Stellen $6.9\cdot10^{-36}$ bei $N=80$* | Experiment |
| 1 | Szekeres' Kennzeichnung (it/br) | Die Abel-Funktion von $e^x-1$ ist die einzige, die total monoton ist (Szekeres 1961) — eine **reelle** Eindeutigkeitsbedingung; Gegenstück zu `br.frage.bohr_mollerup` | Literatur |
| 1 | Shishikura: gestörte Fatou-Koordinaten (br) | Knesers slog = gestörte Fatou-Koordinate; das Werkzeug für `br.frage.implosion` | Literatur |
| 2 | Möbius-Iteration (it) | Bei Möbius-Abbildungen stimmen die regulären Iterationen an beiden Fixpunkten überein — Gegenstück zu „Bummer“ | zwei Zeilen |
| 2 | Grenzfälle (gr) | $\lim_{x\to0^+}{}^nx=1$ für gerade, $0$ für ungerade $n$; ${}^n(-1)=-1$; ${}^{-1}a=0$, ${}^{-2}a=-\infty$; ${}^{-1}1$ beliebig | elementar |
| 2 | Keine Additionsregel (gr) | Es gibt keine Operation $*$ mit ${}^{m+n}x={}^mx*{}^nx$ ($65536$ gegen $256$); Potenz-Identität $({}^ba)^{({}^ca)}=({}^{c+1}a)^{({}^{b-1}a)}$ | elementar |
| 2 | Reihen des endlichen Turms (fp) | Galidakis: $z\uparrow\uparrow m$ in $\ln z$ stimmt bis Ordnung $m$ mit Knoebels Reihe $\sum(n+1)^{n-1}(\ln z)^n/n!$ des unendlichen Turms überein | Koeffizienten vergleichen |
| 2 | Quadratische Näherung (kn) | $1+\frac{2\ln a}{1+\ln a}x-\frac{1-\ln a}{1+\ln a}x^2$ auf $(-1,0]$, $C^1$, nicht $C^2$; ${}^{1/2}2\approx1.45933$ (Kneser $1.45878$) | nachrechnen |
| 2 | Nachbar-Turm $x^{(1/x)^{(1/x)\cdots}}$ (fp) | Konvergiert für $e^{-1/e}\le x\le e^e$, Grenzwert $\ln x/W(\ln x)$, bei $x=e$ die Omega-Konstante | elementar |
| 2 | Integral-Spiegelung (gr) | $\int_0^\infty x^{1/x-2}dx=\int_0^\infty x^{-x}dx=1.995455958$ (Substitution $x\to1/x$) | zwei Zeilen |
| 2 | Ackermann und Wachstum (br) | ${}^n2=A(4,n-3)+3$; Tetration $\approx f_3$ der Fast-Growing-Hierarchy; Anwendung: Fox' Beweis des Graph-Removal-Lemmas mit $2\uparrow\uparrow5h^4\log(1/\varepsilon)$ | Literatur |
| 2 | Ordinalzahlen (br) | $\omega\uparrow\uparrow\omega=\varepsilon_0$; bei $\omega\uparrow\uparrow(\omega+1)$ bricht die naive Fortsetzung | Literatur |
| 2 | Geschichte (Kapitel 0, Textteil) | Goodstein 1947 (Name), Euler 1783, Eisenstein 1844, Barrow 1936, Kneser 1950, Szekeres 1958/61, Tabellen 1969, Asimov–Hickerson–Schroeppel 1991, Walker 1991, Forum ab 2007, Kouznetsov 2009, Paulsen–Cowgill 2017, Paulsen 2019, Nixon 2021 | — |
| 3 | Kommutative und linke Hyperoperationen (br) | $x\triangle_{n+1}y=\exp(\log x\triangle_n\log y)$; linke Klammerung gibt $x^{x^{n-1}}$; balancierte Tetration = Iteration von $x^x$ bzw. $xe^x$ | elementar |
| 3 | Endziffern-Fallen (br) | Die naive Rekursion für die letzten $d$ Ziffern ist falsch, wenn Basis und Modul nicht teilerfremd sind (Gegenbeispiele) | nachrechnen |
| 3 | Erste Ziffern (br) | Führende Ziffern von ${}^ya$ für $y>6$ praktisch unberechenbar (verlangt Ziffern von $\log_{10}a$ an Stelle ${}^{y-2}a$) | — |

## C. Vermutungen und Fragen aus Forum und Wikis (als Vermutung/Frage aufnehmen)

| Vorschlag | Inhalt | Stand |
|---|---|---|
| Basiswechsel nirgends analytisch | siehe A; Gegen-Vermutung zu `bw.verm.streifen` | offen, Forum 2011 |
| TPID 7 | $0<({}^yx)({}^{-y}x)<1$ für $\lvert y\rvert<1$ (regulär unter $\eta$, Kneser darüber); **Stichprobe** $y=\tfrac12$: $0.782$ ($\sqrt2$), $0.821$ ($e$) | offen; leicht flächig prüfbar |
| TPID 10 | Die Folgen der „Eulers“ $(e,e_4,\dots)$ und „Etas“ $(\eta,\eta_4,\dots)$ gehen gegen $4$ bzw. $2$ | offen; abhängig von der Fortsetzung |
| TPID 1 / 3 | Eigenwerte abgeschnittener Carleman-Matrizen gehen gegen $\ln(a)^k$ bzw. $p^k$ | offen; verbindet sich mit `kn.konstr.carleman` |
| TPID 18 | Carleman-Polynom-Tetration gegen Kneser | offen; siehe A |
| TPID 2 | Exponentielle Fakultät $EF(x)=x^{EF(x-1)}$ eindeutig unter totaler Monotonie, dann $EF(0)=\gamma$? | offen, eher Rand |
| TPID 12 / 13 | Newton-Interpolation von $y_n$ (${}^ny_n=n$) bzw. $n^{1/n}$ konvergiert? TPID 13 per Mellin/Ramanujan fast bewiesen | teilweise |
| Rechte Halbebene (Nixon 2017) | Für $1<b<\eta$ bildet die reguläre Tetration $\{\mathrm{Re}z>0\}$ in sich ab ⇒ beschränkte analytische Pentation, Hexation, …; Lücke: Vorzeichen der Taylor-Koeffizienten der Poincaré-Funktion | prüfbar |
| β-Methode (Nixon 2021) | holomorphe Tetration, die **nicht** Kneser ist (nach oben gegen $\infty$ statt $L$); dazu $2\pi i$-periodische Tetrationen | Preprint; zeigt, dass die senkrechte Bedingung nötig ist |
| Ableitungs-Vorzeichen (Fox, Sheldon) | Alle ungeraden Ableitungen von $T_e$ positiv; gerade wechseln genau einmal das Vorzeichen (22. bei $x=0.0078$) | prüfbar mit unseren Ableitungen |
| Fläche der Shell–Thron-Region; Basen mit grösstem Imaginärteil | offene Zahlenfragen (Forum 2022) | in Minuten berechenbar — eigene Werte |
| $\int_{-1}^0T_b=\tfrac12$? | Für welche Basis gilt das bei Kneser (bei linear immer)? | berechenbar |
| Offene Rationalitätsfragen (ar) | Ist ${}^4\pi$ oder ${}^5e$ eine ganze Zahl? Ist ${}^nq$ ($q$ rational, nicht ganz) je rational? Ist die positive Lösung von ${}^nx=2$ ($n>3$) rational? | offen (Wikipedia) |

## D. Bewusst nicht aufnehmen

- Notations-Zoo (Bowers, BEAF, Hyper-E, Plus-Notation …) und Grosszahlen-Tabellen — kein mathematischer Inhalt über Tetration.
- Spekulative Foren-Fäden ohne prüfbaren Kern (Zeration, „fake function theory“, Collatz/Primzahl-Exkurse, „Tetrion space“, Neuronale Netze gegen Kneser).
- TPID 16/17 (zu vage formuliert), „Tetration ist Pentation“-Spekulationen.
- Physik-Analogien (Turbulenz, Hopf-Faserung).

## Vorschlag zum Vorgehen

1. **Runde 1 (Ehrlichkeit): erledigt 08.10.2026.** Alle Punkte aus A eingetragen: 15 neue Quellen (Forum, Wikipedia, Mizugadro, Euler,
   Eisenstein, Barrow, Macintyre, Knoebel); Neuheit „bekannt“ für `kl.beob.zwei_abel`, `bw.beob.m2`, `kl.satz.eta_koordinate`; Vorarbeiten in
   `bw.def.F_bd`, `bw.satz.phasengesetz`, `bw.prop.hub`; Hooshmand im Wortlaut; Walker als Basiswechsel-Grenzwert; neue Knoten
   `bw.verm.nirgends_analytisch` (Gegen-Vermutung) und `kl.frage.eta_analytisch`; `bw.verm.streifen` als umstritten markiert.
   Nebenbei gefunden und behoben: Der Quellenschlüssel `justus_p2` stand doppelt (Paper II und Periode-2-Notiz), YAML hatte den ersten
   stillschweigend überschrieben — jetzt `justus_periode2`, und der Lader meldet doppelte Schlüssel als Fehler.
2. **Runde 2 (Kern): erledigt 08.10.2026.** 13 neue Aussagen mit Theorie-Block, vier Experimente (`exp.superwurzel` 10/10, `exp.tpid6` 4/4,
   `exp.qbinomial` 10/10, `exp.cheta` 13/13), sechs neue Quellen (Szekeres 1961, Hyperoperations-Wiki, Shishikura, Lavaurs, Reshetnikov MO 259278
   im Original, Fox „cheta“), neuer Strang `gr.umkehrungen`. Befunde: TPID-6-Rate $n^2(x_n-\eta)\to2\pi^2\eta/e$ auf 5 Stellen
   (`fp.beob.tpid6_rate`, vermutlich neu); die q-Binomial-Reihe konvergiert doch — die „wachsenden Glieder“ in Tabelle B waren
   Auslöschung bei zu wenig Stellen (nötig $\approx\tfrac{n^2}2\log_{10}(1/q)$); Rate $q^{\mathrm{Re}z+2}$ gemessen und erklärt;
   Reshetnikovs Vermutung 3 ist gleichwertig zu „alle Koenigs-Koeffizienten negativ“, für vier Basen bis $c_{60}$ bestätigt (`kl.verm.vollmonoton`);
   Szekeres' total monotone Abel-Funktion bis zur 8. Ableitung nachgerechnet. Die Ansicht erkennt jetzt Listen direkt nach einem
   Einleitungssatz (71 Stellen im Lehrbuch waren bisher Fliesstext).
3. **Runde 3 (Höhere Stufen): erledigt 08.10.2026.** Neuer Strang `kn.hoehere_stufen` (Kapitel 6) mit drei Aussagen, zwei Experimenten
   (`exp.pentation` 9/9, `exp.tetra_euler` 10/10) und drei Forumsquellen (Fox 2007 tid 103; tid 262 mit Levensteins 32 Stellen von 2010;
   Levenstein 2020 tid 1263). Pentation zur Basis $e$ am Fixpunkt $L=-1.85035452902718141848\ldots$ ($\lambda=6.46067129568183939\ldots$):
   Fox' Werte (Matrix-slog, 2007) stimmen mit Kneser auf 22 Stellen; $e\uparrow\uparrow\uparrow\tfrac12=1.6323247404\ldots$; Basis 2 trifft
   Levensteins $2\uparrow\uparrow\uparrow1.5\approx2.7029$. $\eta_4$ und $e_4$ auf alle 32 Stellen von Levenstein bestätigt (Abstand $3\cdot10^{-32}$);
   die drei Kennzeichnungen (Fixpunktzahl 3 → 1, Maximum der Selbst-Tetrawurzel, Konvergenz von $b\uparrow\uparrow\uparrow\infty$) nachgerechnet;
   neu dabei $1.6\uparrow\uparrow\uparrow\infty=2.26438440348\ldots$ und $\sqrt2\uparrow\uparrow\uparrow\infty=1.54921732984\ldots$ Hinweis: Der Wert
   $3.08855\ldots$ aus der Matrix-Rechnung 2009 ist ungenau, nicht falsch zitiert — er hängt an der Näherung.
4. **Runde 4 (Vermutungen und Rand): erledigt 08.10.2026.** 12 Aussagen, vier Experimente (`exp.vermutungen` 9/9, `exp.shell_thron` 6/6,
   `exp.ableitungen_vorzeichen` 6/6, `exp.stufe5` 7/7), acht Quellen, Geschichte als Abschnitt in Kapitel 0, Tabellen im Lehrbuch. Befunde:
   - **TPID 7** hält regulär und für Kneser bis $b_7=10.2210634157814\ldots$ (dort $T''(0)=T'(0)^2$); darüber **Gegenbeispiel** ($b=100$: $1.316$).
     Ursache: Log-Konkavität von $T$.
   - **Halbebene (Nixon 2017):** Aus „alle Koenigs-Koeffizienten negativ“ (Reshetnikov 3) folgen Konkavität, TPID 7 und $|T(z)-L|\le L-1$ auf
     $\mathrm{Re}z\ge0$ (`kl.satz.folgen_vollmonoton`, bedingt); gemessen für vier Basen.
   - **TPID 10:** $\eta_5=1.7373559251516875731\ldots$, $e_5=3.3014072785233674934\ldots$ (vermutlich erstmals berechnet); beide Folgen wachsen, TPID 11 gilt.
   - **TPID 1:** Carleman-Eigenwerte → $\lambda^k$ (Basis 1.3: $3\cdot10^{-16}$ bei $N=32$).
   - **Ableitungen (mike3/Levenstein 2014):** ungerade bis 25 positiv; gerade wechseln je einmal, $x_{22}=0.0077902$ (der Wert in Tabelle C stimmt).
   - **Shell–Thron:** Fläche $8.0509332809395\ldots$ und Extrempunkte — schon von Robbins 2009 (tid 263) auf 30 Stellen berechnet, hier unabhängig bestätigt.
   - $\int_{-1}^0T_b=\tfrac12$ bei $b=2.6937028415\ldots$; ${}^2q$ irrational, ${}^3q$ transzendent für rationale nicht ganze $q$; β-Methode als offene Frage
     (Nixon zog die Behauptung „neue Tetration“ selbst zurück).
   Nicht aufgenommen: TPID 2, 3, 12, 13 (Rand bzw. ohne eigene Prüfung) — TPID 2, 12, 13 in Runde 5 nachgeholt.
5. **Runde 5 (Vollabgleich): erledigt 10.10.2026.** Diesmal das ganze Forum: alle 1265 Themen der Fachforen als Druckfassung gespiegelt (ein Abruf pro
   Sekunde, nur lokal ausgewertet), alle Titel gesichtet, die rund 620 Themen über 4 kB über Kopf und ersten Beitrag, rund 70 genau gelesen.
   **32 neue Aussagen, sieben Experimente** (`exp.ln_basis`, `exp.turm_pde`, `exp.forum_formeln`, `exp.halbe_exp`, `exp.zwei_sinh`, `exp.forum_rand`,
   `exp.forum_iteration`), alle bestanden:
   - *Ableitungen und Basis:* ln b = T′(0)/T′(−1) = slog′(0)/slog′(1) lösungsunabhängig (hyper4geek 2007, Robbins 2009); slog′ zwischen allen Gliedern
     der Logarithmen-Skala; die Basisableitung als lineare Stufengleichung samt Summenformel bei ganzer Höhe; Kobis PDE für Türme (tid 417, 2010) —
     **eigener Befund:** Kneser erfüllt sie bei gebrochener Höhe nur bis auf ein periodisches Wackeln (Spanne 5.7·10⁻⁴ bei Basis 2), die infinitesimale
     Fassung der Basiswechsel-Welle.
   - *Halbe Exponentialfunktion* (Knesers Ausgangsfrage, tid 228/544/1043/1179): Werte, Nullstelle, Singularitäten, Wachstum zwischen Polynomen und e^{x^ε}.
   - *Andere Konstruktionen:* 2sinh-Methode (gegen Kneser nachgerechnet: 3·10⁻⁴ daneben), Kontinuumssummen (Ansus, mike3), Kouznetsovs Cauchy-Integral,
     quadratische Näherung, der allgemeine Grenzwert-Trick hinter Basiswechsel, Walker, 2sinh und β.
   - *Formeln:* rationale Potenzsummen der Fixpunkte von e^z (Fox 2007), Newton-Reihe der ganzen Iterierten = reguläre Tetration (bewiesen), x↑↑n und der
     unendliche Turm bis (ln x)^n gleich, i^i^i^…, Nachbar-Turm, Grenzfälle bei 0, keine Additionsregel.
   - *Iteration und Eindeutigkeit:* elementare Superfunktionen, Möbius-Modell für den Übergang bei η, Carlson-Eindeutigkeit unter η (mike3 2010),
     TPID 4 (eigener Zweizeilen-Beweis), TPID 9 (Trappmann 2010), TPID 2, 5, 12/13 als Fragen; „Grand Unity Conjecture“ als Quelle bei der Carleman-Vermutung.
   - *Rand:* Ackermann/schnell wachsende Hierarchie, Ordinalzahlen (ε₀), linke/ausgewogene/Bennett-Operationen, exponentielle Fakultät (Liouville-Zahl),
     Vergleich gebrochener Iterierter verschiedener Basen (Taranovsky/Levenstein).
   **Korrektur nebenbei:** b₇ aus TPID 7 ist 10.2210634157814…, nicht 10.2210634581… (zu grobe Abtastung der Ableitungen bei grossen Basen); e34 korrigiert.
   **Bewusst nicht aufgenommen** (Grund): Zeration/Deltation, Nept/Nopt, Notationen, Salat-Zahlen (kein mathematischer Inhalt über Tetration);
   tommys Spekulationen ohne prüfbaren Kern (Fake-Funktionen/TPID 17, TPID 14–16, Gauss-Methode im Detail, Spiralzahlen, Zahlentheorie-Exkurse);
   Helms' alternierende Iterationsreihen (Stand unklar, Gegenbeispiel von G. A. Edgar); Nixons Bruchrechnungs-Ansatz und Halboperationen (unvollständig;
   als Frage im Bennett-Knoten); p-adische Fragen (über Ripà hinaus nichts Gesichertes); Superfunktionen anderer Funktionen (Fakultät, zex); Levensteins
   Konstante zwischen den beiden η-Superfunktionen (ungeprüft); natürliche Grenze der regulären Tetration in der Basis (offen, unklar formuliert); TPID 3.

Jede neue Aussage bekommt wie bisher einen Theorie-Block, Quellen mit Prüfvermerk und — wo möglich — ein Experiment.

## Quellen

- https://en.wikipedia.org/wiki/Tetration · https://en.wikipedia.org/wiki/Iterated_logarithm · https://en.wikipedia.org/wiki/Iterated_function
- https://de.wikipedia.org/wiki/Potenzturm · https://de.wikipedia.org/wiki/Hyper-Operator · https://ru.wikipedia.org/wiki/Тетрация
- https://googology.fandom.com/wiki/Tetration · https://mathworld.wolfram.com/PowerTower.html
- https://web.archive.org/web/20241206205947/https://mizugadro.mydns.jp/t/index.php/Tetration
- https://tetrationforum.org/ — FAQ (tid 184), offene Probleme (tid 162), Forums-Wiki https://tetrationforum.org/hyperops_wiki/
- Forum-Themen: tid 14, 26, 69, 103, 108, 236, 260, 316, 333, 360, 366, 372, 476, 578, 648, 649, 825, 1143, 1168, 1182, 1188, 1334, 1354, 1570, 1668, 1824
- Runde 5: alle 1265 Themen (Spiegel nur lokal, nicht im Repo); zitiert u. a. tid 10, 12, 47, 91, 94, 120, 124, 145, 150, 206, 214, 275, 337, 341, 343, 353, 361, 417, 424, 437, 475, 496, 544, 719, 789, 809, 819, 1043, 1053, 1144, 1179, 1189, 1290, 1559, 1608; arXiv:1008.1409 (Trappmann)

# Literaturkarte fp — Fixpunkte, unendlicher Turm, Parameterraum der Exponentialfamilie (L1, Stand 12.10.2026)

Prüfvermerk wie in KARTE_it.md (O Original · R Rezension · A Abstract · Z Zitat/Wikipedia/Forum · G Gedächtnis).

Bestand vorher: 41 Aussagen (fp.*) in 03_ (Fixpunkte, Turm) und 09_ (Parameterraum, Periode 2, Mittelwerte, gleiches L).
Literatur: Euler, Eisenstein, Barrow, Thron, Shell, Macintyre, Knoebel, Baker–Rippon (Fibonacci 1985), Corless u. a., Misiurewicz, Devaney,
Schleicher. Die Periode-2- und Siegel-Knoten waren rein eigen.

## 1. Unendlicher Turm und Lambert W

| Jahr | Autor(en) | Werk, Ort — Schlüssel | Hauptergebnis | Methode | Beweisstand | Prüf. | Bei uns |
|---|---|---|---|---|---|---|---|
| 1758 | J. H. Lambert | Observationes variae…, Acta Helvetica 3 — `lambert1758` | Reihe der trinomischen Gleichung (Ursprung von W) | Reihen | formal | G | — |
| 1778 | L. Euler | De formulis exponentialibus replicatis (E489), Acta Petropol. 1777:I — `euler1778` | Konvergenz von $r, r^\alpha, r^{r^\alpha},\dots$ (Frage Condorcets); Grenze $e^{1/e}$ | Fixpunktbetrachtung | für die Grenze | Z (Euler-Archiv) | fp.satz.turm_reell (jetzt zitiert) |
| 1783 | L. Euler | De serie Lambertina (E532) — `euler1783` | Lambertsche Reihe, Vorform von W | Reihen | | Z | fp.satz.turm_reell, fp.satz.lambert |
| 1844 | G. Eisenstein | Entwicklung von $\alpha^{\alpha^{\cdots}}$, J. reine angew. Math. 28, 49–52 — `eisenstein1844` | Reihe des Grenzwerts in $\ln\alpha$ | Reihen | | Z | fp.satz.turm_reell, fp.satz.turm_reihe_ln |
| 1936 | D. F. Barrow | Infinite exponentials, AMM 43, 150–160 — `barrow1936` | Konvergenz genau für $e^{-e}\le a\le e^{1/e}$, darunter Pendeln zwischen zwei Werten; veränderliche Basen | reell | bewiesen | R (JFM) | fp.satz.turm_reell, fp.satz.p2_form (jetzt) |
| 1945 | J. Ginsburg | Iterated exponentials, Scripta Math. 11 — `ginsburg1945` | Übersicht | | | G | — |
| 1957 | W. J. Thron | Convergence of infinite exponentials with complex elements, PAMS 8 — `thron1957` | Konvergenz, wenn alle $|\ln a_n|\le1/e$ | Abschätzung | bewiesen | Z (Baker–Rippon 1989) | fp.satz.shell_thron |
| 1959 | E. M. Wright | Solution of $ze^z=a$, Bull. AMS 65 — `wright1959` | Lösungen von $ze^z=a$ (Zweige von W) | | bewiesen | G | — |
| 1962 | D. L. Shell | On the convergence of infinite exponentials, PAMS 13 — `shell1962` | Konvergenzgebiet $b=\exp(te^{-t})$, $|t|<1$ (Shell–Thron) | Multiplikator | bewiesen | Z | fp.satz.shell_thron |
| 1966 | A. J. Macintyre | Convergence of $i^{i^{\cdots}}$, PAMS 17 — `macintyre1966` | Turm zur Basis $i$ konvergiert | | bewiesen | Z | fp.beob.i_turm |
| 1980/81 | Creutz, Sternheimer | Fibonacci Quart. 18, 19 — `creutz_sternheimer1980` | Konvergenz iterierter Potenzen, Rechnungen | | | G | — |
| 1981 | R. A. Knoebel | Exponentials reiterated, AMM 88 — `knoebel1981` | grosse Übersicht mit Geschichte und ~60 Literaturangaben; Reihe des Turms; $x^y=y^x$ | Übersicht | | Z | mehrere fp.* |
| 1983 | P. J. Rippon | Infinite exponentials, Math. Gaz. 67 — `rippon1983` | elementare Darstellung | | | G | — |
| 1983–89 | I. N. Baker, P. J. Rippon | Ann. Acad. Sci. Fenn. 8 (1983); Fibonacci Quart. 23 (1985) — `baker_rippon1985`; AMM 92 (1985) — `baker_rippon1985b`; Math. Proc. Camb. 105 (1989) — `baker_rippon1989a`; Complex Var. 12 (1989) — `baker_rippon1989b` | (i) Turm konvergiert gegen einen Grenzwert $\ne$ Glieder nur im Abschluss der Kardioide, innen stets; (ii) zwei abwechselnde Basen: $\log a=se^{-t}$, $\log b=te^{-s}$, Konvergenz ⇔ $|st|\le1$; (iii) periodische Basen allgemein; (iv) veränderliche Basen, grössere Mengen als Barrow/Thron | Fixpunkt + Multiplikator, Normalfamilien | bewiesen | R | fp.satz.shell_thron, fp.satz.p2_form, fp.satz.p2_kein_kreis (jetzt) |
| 1986 | J. M. de Villiers, P. N. Robinson | The interval of convergence and limiting functions of a hyperpower sequence, AMM 93 — `devilliers_robinson1986` | Konvergenzintervall mit schärferen Abschätzungen; Grenzfunktionen reell analytisch, Randwerte | reell | bewiesen | R | fp.satz.turm_reell (jetzt) |
| 1989 | J. MacDonnell | IJMEST 20 — `macdonnell1989` | gerade/ungerade Teilfolgen auf $(0,1)$ | | | R | — |
| 1994/96 | Elstrodt; Länger; Meyerson; Ash | Iterierte Potenzen (Math. Semesterber. 41); elementarer Beweis (Elem. Math. 51) — `laenger1996` | elementare Konvergenzbeweise, Grenzwerte | | bewiesen | R (Länger) | fp.satz.turm_reell (Länger jetzt) |
| 1995 | G. Bachman | Convergence of infinite exponentials, Pacific J. 169 — `bachman1995` | zwei scharfe Konvergenztests für veränderliche Basen | | bewiesen | R | — |
| 1996 | Corless, Gonnet, Hare, Jeffrey, Knuth | On the Lambert W function, Adv. Comput. Math. 5 — `corless1996` | Standardreferenz W; Turm $=-W(-\ln b)/\ln b$ | | bewiesen | G | fp.satz.lambert |
| 1998 | A. P. Bulanov | Regularity of infinite exponentials, Izv. Math. 62 — `bulanov1998` | Holomorphie des Grenzwerts als Funktion der Glieder | | bewiesen | R | — |
| 2003 | S. R. Finch | Mathematical Constants §6.11 (Bib?) | Konstanten iterierter Exponentiale | | | G | — |
| 2004 | Y. Cho, Y.-O. Kim | Analytic properties of the limits of the even and odd hyperpower sequences, Bull. KMS 41 — `cho_kim2004` | $h_e,h_o$ mit $x^{h_e}=h_o$, $x^{h_o}=h_e$; Asymptotik an $e^{-e}$ und an 0 | reelle Analysis | bewiesen | R | fp.satz.p2_form, fp.satz.p2_raten (jetzt, als möglicher Vorläufer) |
| 2004 | J. Anderson | Iterated exponentials, AMM 111 — `anderson2004` | Übersicht, überraschende Antworten für $a>0$ | | | G | — |
| 2004–07 | I. N. Galidakis | Complex Var. 49 (2004), 50 (2005), 52 (2007) — `galidakis2004`, `galidakis2005` | W und Turm; Lösung der 1. und 2. Hilfsgleichung (Fixpunkt, 2-Zyklus) mit W; periodische Punkte; Hyper-Lambert-Funktionen | W-Zweige | bewiesen | R | fp.satz.p2_form (jetzt) |
| 2012 | Gao, Han, Schilling | On the rate of convergence of iterated exponentials, J. Appl. Math. Comput. 39 — `gao_han_schilling2012` | Konvergenzrate des Turms (Inhalt nicht zugänglich) | | ? | — | fp.satz.p2_raten (jetzt, Inhalt offen) |
| 2019 | L. Moroni | The strange properties of the infinite power tower, arXiv — `moroni2019` | didaktische Übersicht | | | A | — |
| Forum | Robbins 2009; tid 94, 120, 1053 | `forum_robbins2009_shellthron` u. a. | Fläche Shell–Thron; Summen über Fixpunkte; $i^{i^{\cdots}}$ | | | O/Z | fp.satz.shell_thron_flaeche, fp.satz.fixpunkt_summen |

## 2. Parameterraum der Exponentialfamilie ($b^z\simeq e^w+\kappa$, $\kappa=\ln\ln b$)

| Jahr | Autor(en) | Werk, Ort — Schlüssel | Hauptergebnis | Beweisstand | Prüf. | Bei uns |
|---|---|---|---|---|---|---|
| 1926 | P. Fatou | Acta Math. 47 — `fatou1926` | Iteration ganzer transzendenter Funktionen; Frage $J(e^z)=\mathbb C$? | bewiesen | R | — |
| 1981 | M. Misiurewicz | On iterates of $e^z$, ETDS 1 — `misiurewicz1981` | $J(e^z)=\mathbb C$ | bewiesen | G | fp.satz.abstossend |
| 1984 | R. L. Devaney; Devaney–Krych | Bull. AMS 11; Dynamics of exp(z), ETDS 4 — `devaney1984`, `devaney_krych1984` | $J(\lambda e^z)=\mathbb C$ für $\lambda>1/e$, Explosion bei $1/e$; Symbolik | bewiesen | R | fp.satz.abstossend, fp.beob.regionen |
| 1984 | I. N. Baker, P. J. Rippon | Iteration of exponential functions, Ann. Acad. Sci. Fenn. 9 — `baker_rippon1984` | keine wandernden Gebiete für $e^{az}$; $D_1$, $D_2$ (anziehende 1- und 2-Zyklen) zusammenhängend; $J=\mathbb C$ für reelle $a>1/e$ und dicht auf $\partial D_1$ | bewiesen | R | fp.satz.p2_form (jetzt); eigener Knoten fehlt |
| 1985 | M. R. Herman | Are there critical points on the boundaries of singular domains?, CMP 99 — `herman1985` | diophantische Drehzahl ⇒ Siegel-Scheibe der Exponentialabbildung unbeschränkt | bewiesen | Z (Rempe 2004) | fp.satz.siegel_unbeschraenkt (neu) |
| 1986/87 | M. Rees; M. Lyubich | Math. Z. 191; Sib. Math. J. 28 — `rees1986`, `lyubich1987` | messbare Dynamik von exp: nicht rekurrent, Bahnen fast aller Punkte häufen sich am Orbit von 0 | bewiesen | R | — |
| 1987 | C. McMullen | Trans. AMS 300 — `mcmullen1987` | $J(\lambda e^z)$ hat Hausdorff-Dimension 2; bei anziehendem Zyklus Fläche 0 | bewiesen | R | — |
| 1987 | Devaney, Goldberg | Duke 55 — `devaney_goldberg1987` | Uniformisierung anziehender Becken | bewiesen | R | — |
| 1992 | Eremenko, Lyubich | Ann. Inst. Fourier 42 — `eremenko_lyubich1992` | Klasse $\mathcal S$: keine wandernden Gebiete, Fluchtmenge | bewiesen | G | fp.satz.p2_kein_kreis (nur Wort) |
| 1999 | Bodelón, Devaney, Hayes, Roberts, Goldberg, Hubbard | IJBC 9 — `bodelon1999` | Haare in Dynamik- und Parameterebene | bewiesen | R | — |
| 2002 | Devaney, Fagella, Jarque | Fund. Math. 174 — `devaney_fagella_jarque2002` | jede S-Knetfolge wird von einer hyperbolischen Komponente realisiert | bewiesen | R | — |
| 2003–09 | D. Schleicher; Zimmer; Rempe; Förster | `schleicher2004`, `schleicher_zimmer2003`, `rempe_schleicher2003` (= Invent. Math. 175, 2009), `rempe_schleicher2008`, `foerster_rempe_schleicher2008`, `foerster_schleicher2009` | hyperbolische Komponenten: Rand zusammenhängend, konforme Uniformisierung auf Halbebene; Bifurkationslocus zusammenhängend; Parameterstrahlen; Kombinatorik der Kind-Komponenten | bewiesen | A | fp.satz.exp_normalform, fp.satz.p2_kein_kreis; fp.verm.knospen (Kind-Komponenten!) |
| 2004 | L. Rempe | On a question of Herman, Baker and Rippon concerning Siegel disks, Bull. LMS 36 — `rempe2004` | unbeschränkte Siegel-Scheibe ⇒ singulärer Wert auf dem Rand | bewiesen | **O** (Einleitung und Satz) | fp.satz.siegel_unbeschraenkt (neu) |
| 2001 | L. Geyer | Trans. AMS 353 — `geyer2001` | Brjuno notwendig und hinreichend für Linearisierung von $\lambda ze^z$ (andere Familie) | bewiesen | R | fp.satz.neutral |

## 3. Dieselbe Frage, verschiedene Zugänge

**F1 — Wann konvergiert der Turm?** reell monoton (Euler, Barrow, de Villiers–Robinson, Länger: elementar); Fixpunkt + Multiplikator (Thron, Shell, Baker–Rippon 1985 AMM: innen ja, aussen nein); geschlossene Form über W (Eisenstein, Corless u. a.); komplexe Dynamik (Baker–Rippon 1984: $D_1$ als hyperbolische Komponente); veränderliche/periodische Basen (Thron, Baker–Rippon 1983/89, Bachman).

**F2 — Was passiert unter $e^{-e}$?** Barrow (Pendeln), Cho–Kim (Analytik der beiden Grenzfunktionen), Baker–Rippon 1985/89 (Parametrisierung durch $|st|\le1$), Galidakis 2005 (W-Lösung des 2-Zyklus), Baker–Rippon 1984 ($D_2$ zusammenhängend). Eigene Fassung: fp.satz.p2_form ($G(p)$-Form), fp.satz.rechteckwelle (Tetration statt Turm — dafür keine Literatur gefunden).

**F3 — Was geschieht an irrationalen Randpunkten?** Siegel/Brjuno (Scheibe existiert), Cremer (generisch keine), Herman (unbeschränkt), Rempe (Turm auf dem Rand). Eigener Mittelwertsatz fp.satz.siegel_mittel war auf eine falsche Annahme gebaut → jetzt Skizze.

## 4. Befunde zur Neuheit (12.10.2026)

- **fp.ann.siegel_rand widerlegt** (Herman 1985): für diophantisches α — auch den goldenen Schnitt, an dem gerechnet wurde — ist die Siegel-Scheibe unbeschränkt. Neuer Literaturknoten fp.satz.siegel_unbeschraenkt; fp.satz.siegel_mittel von „bedingt“ auf „skizze“ (Lücken benannt).
- **fp.satz.p2_form**: im Kern bekannt (Barrow 1936 reell; Baker–Rippon 1985/89 Parametrisierung; Baker–Rippon 1984 $D_2$; Galidakis 2005) → Neuheit „bekannt“.
- **fp.satz.p2_raten**: Standardrechnung am parabolischen Punkt; mögliche Vorläufer der Konstanten Cho–Kim 2004, Gao–Han–Schilling 2012 → „folklore“.
- **fp.satz.turm_reell**: Euler-Zuschreibung präzisiert — die Konvergenzfrage steht in E489 (1778); E532 (1783) ist die Lambert-Reihe. MathWorld schreibt das ganze Intervall „Euler 1783“ zu; welche der beiden Arbeiten die untere Grenze $e^{-e}$ enthält, ist ungeprüft (Mainzer Übersetzung von E489 lesen).
- **fp.beob.tpid6_rate** (vermutlich neu): Methode klassisch (Engpass, Pomeau–Manneville; parabolische Rate), Anwendung auf TPID 6 nicht gefunden; Gao–Han–Schilling unbekannt → bleibt, Lücke 3.

## 5. Lücken (Priorität)

1. erledigt 12.10.: fp.verm.knospen jetzt Satz (Rempe–Schleicher 2008, Prop. 6.2/Thm. 6.8(c)) — fp.verm.knospen ist vermutlich ein Satz: Kind-Komponenten an Randpunkten mit Multiplikator $e^{2\pi ip/q}$ (Rempe–Schleicher 2008/2009, Schleicher 2003 „Attracting dynamics“) — Satzstelle im Volltext finden; anziehender Zyklus zieht den singulären Wert an (Fatou) ⇒ Turm läuft hinein.
2. erledigt (fp.satz.baker_rippon_parameter, fp.satz.baker_rippon_kardioide, fp.satz.baker_rippon_periodisch) — Baker–Rippon 1984 als Knoten ($D_1$, $D_2$ zusammenhängend, keine wandernden Gebiete); Baker–Rippon 1985 (AMM) als eigener Satz „nur im Abschluss der Kardioide“.
3. **hoch** — Volltexte Cho–Kim 2004 und Gao–Han–Schilling 2012 (Raten, Konstanten) für fp.satz.p2_raten und fp.beob.tpid6_rate.
4. erledigt (fp.satz.grenzfunktionen, fp.satz.thron_variabel) — de Villiers–Robinson 1986 (Grenzfunktionen analytisch), Thron 1957 / Baker–Rippon 1989 (veränderliche Basen), Bachman 1995 als Knoten.
5. erledigt (fp.satz.messbare_dynamik) — Messbare Dynamik und Dimension (Rees, Lyubich, McMullen) als Hintergrund zu fp.satz.abstossend / fp.beob.regionen.
6. **mittel** — Knoebel 1981: Literaturliste (~60 Titel) einsehen (JSTOR, nicht frei) — wichtigste Bibliografie des Turms.
7. **gering** — Lambert 1758, Wright 1959, Ginsburg 1945, Creutz–Sternheimer, Anderson 2004, Galidakis' Hyper-Lambert-Funktionen, Moroni 2019.

## 6. Stand und Restlücken (Abschluss L1, 12.10.2026)

- **Stand:** 41 → 48 Aussagen (fp.*). Neu (`daten/aussagen/lit_fp_1.yaml` und 09_): fp.satz.siegel_unbeschraenkt, fp.satz.baker_rippon_kardioide,
  fp.satz.baker_rippon_periodisch, fp.satz.baker_rippon_parameter, fp.satz.thron_variabel, fp.satz.grenzfunktionen, fp.satz.messbare_dynamik.
  fp.verm.knospen ist jetzt Satz (Rempe–Schleicher 2008); fp.ann.siegel_rand widerlegt (Herman 1985); fp.satz.p2_form „bekannt“, fp.satz.p2_raten „folklore“,
  fp.satz.rechteckwelle nach dokumentierter Recherche „vermutlich neu“ (Turm-Literatur und Forum tid 352 kennen die Rechteckwelle nicht).
- **Restlücken:** Volltexte Cho–Kim 2004 (nur PDF) und Gao–Han–Schilling 2012 (Bezahlschranke) — entscheiden über die Neuheit der Konstanten in
  fp.satz.p2_raten und fp.beob.tpid6_rate; Knoebel 1981 (JSTOR) mit seiner Literaturliste; Baker–Rippon 1983 (keine Rezension); welche Euler-Arbeit die
  untere Grenze $e^{-e}$ enthält (E489 in der Mainzer Übersetzung lesen); Cremer-Punkte für $b^z$ und Notwendigkeit der Brjuno-Bedingung in der Exponentialfamilie;
  Wortlaut von Rempe–Schleicher Prop. 6.2 vollständig lesen.

**Nachtrag 12.10. 18:15 (freie PDFs gelesen):** Euler E489 (Mainzer Übersetzung, §§38–39) hat schon **beide** Grenzen und das Pendeln zwischen zwei Werten
unter $1/e^e$ — die Zuschreibung „Barrow 1936“ für das Pendeln ist damit ersetzt (fp.satz.turm_reell, fp.satz.p2_form). Cho–Kim 2004 im Original:
Thm. 2.1 entwickelt $h_e,h_o$ in $\sqrt{e^{-e}-x}$ mit Leitkoeffizient $\sqrt{6/e}$ — unsere Lücke $2\sqrt6e^{-3/2}\sqrt t$ in fp.satz.p2_raten ist damit
bekannt; die Rate $\tfrac{\sqrt6}en^{-1/2}$ der Folge steht dort nicht. Knoebel 1981 und Gao–Han–Schilling 2012 sind nicht frei verfügbar (Restlücke).

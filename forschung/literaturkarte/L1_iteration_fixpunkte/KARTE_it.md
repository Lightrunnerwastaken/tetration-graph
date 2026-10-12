# Literaturkarte it — lokale Iterationstheorie (L1, Stand 12.10.2026)

**Prüfvermerk** (Spalte „Prüf.“): Bibliografie jeweils über zbMATH Open (API), Crossref oder Verlagsseite geprüft, ausser „Bib?“.
Inhalt: **O** Original gelesen · **R** Rezension (zbMATH/JFM) · **A** Abstract · **Z** zitierende Arbeit / Wikipedia / Forum · **G** aus dem Gedächtnis (ungeprüft).
„Bei uns“: Knoten, die das Werk schon tragen; **—** = Lücke. Schlüssel = Eintrag in `daten/quellen/lit_L1.yaml` bzw. `daten/quellen.yaml`.

Bestand vorher: 12 Aussagen (it.*), davon 8 klassisch/bekannt, 4 mit eigenem Anteil; Quellen fast nur Milnor, Kuczma, Szekeres, ein unscharfer Baker-Schlüssel.

## 1. Werke (chronologisch)

| Jahr | Autor(en) | Werk, Ort — Schlüssel | Hauptergebnis | Methode | Beweisstand | Prüf. | Bei uns |
|---|---|---|---|---|---|---|---|
| 1815 | C. Babbage | Essay towards the calculus of functions, Phil. Trans. 105 — `babbage1815` | erste systematische Behandlung von $f(f(x))=x$ und funktionalen Wurzeln; Konjugation $\varphi^{-1}g\varphi$ erzeugt neue Lösungen | Konjugation | vollständig für Beispiele | Z | — |
| ~1824/1881 | N. H. Abel | Détermination d'une fonction…, Œuvres II — `abel1881` | Abel-Gleichung $\alpha(f(x))=\alpha(x)+1$: Iteration wird Verschiebung | Funktionalgleichung | formal | G (Titel geprüft) | fg.def.abel (nur Kuczma) |
| 1871 | E. Schröder | Ueber iterirte Functionen, Math. Ann. 3 — `schroeder1871` | Schröder-Gleichung $\sigma\circ f=\lambda\sigma$; explizite Iterierte gebrochen linearer und trigonometrisch/elliptisch verkleideter Abbildungen; logistische Abbildung $4x(1-x)$ (S. 306) | Konjugation auf linear | vollständig für die Beispiele | R, Z | it.satz.moebius_iteration, it.satz.moebius_familie, it.satz.elementare_superfunktionen (jetzt zitiert) |
| 1884/85 | G. Koenigs | Recherches sur les intégrales…, Ann. ÉNS — `koenigs1884`, `koenigs1885` | Linearisierung am hyperbolischen Fixpunkt, $\sigma=\lim\lambda^{-n}(f^n-p)$; Abel-Lösung mit log. Singularität eindeutig | Grenzwert, Normalität | bewiesen | R (1885) | it.satz.koenigs |
| 1897 | L. Leau | Étude sur les équations fonctionnelles, Ann. Toulouse 11 — `leau1897` | Multiplikator Einheitswurzel / 1: Blütenblätter, formale Normalformen | Potenzreihen, Sektoren | bewiesen | R | it.satz.leau_fatou (zitiert nur Milnor) |
| 1904 | L. E. Böttcher | Izv. Kazan 14 — `boettcher1904` | superanziehender Fixpunkt: Koordinate mit $\phi(f)=\phi^k$ | Grenzwert $\lim (f^n)^{1/k^n}$ | bewiesen | G (Bib? über Wikipedia) | — |
| 1918 | G. Julia | Mémoire…, J. Math. Pures Appl. — `julia1918` | globale Theorie rationaler Iteration (Julia-Menge); Julia-Gleichung für die Erzeugende eines Flusses | Normalfamilien | bewiesen | R | it.lem.julia (nur Kuczma) |
| 1919–20 | P. Fatou | Sur les équations fonctionnelles, Bull. SMF — `fatou1919` | Fatou-Koordinaten, Blume (mit Leau), Fatou/Julia-Mengen | Normalfamilien | bewiesen | R | it.satz.leau_fatou |
| 1923–24 | P. Fatou | Itération analytique et substitutions permutables — `fatou1924` | vertauschbare Abbildungen und Koenigs-Iteration | Koenigs-Koordinaten | bewiesen | R | it.satz.vertauschbar (nur Kuczma) |
| 1927–38 | H. Cremer | Zum Zentrumproblem, Math. Ann. 98; Häufigkeit der Nichtzentren, Math. Ann. 115 — `cremer1928`, `cremer1938` | es gibt irrational indifferente Fixpunkte ohne Linearisierung (Cremer-Punkte); generisch | Kleine Nenner, Kategorie | bewiesen | R | fp.satz.neutral (nur als Wort) |
| 1927–28 | P. Lévy | Sur l'itération de la fonction exponentielle, C. R. 184; Fonctions à croissance régulière…, Ann. Mat. 5 — `levy1927`, `levy1928` | erste gebrochene Iteration von $e^x$: Rückführung auf $e^x-1$ am parabolischen Punkt 0; „reguläres Wachstum“ als Auswahlprinzip | Konjugation ∞ ↔ parabolischer Punkt | bewiesen (Existenz), Eindeutigkeit über Wachstum | R (Noten) | — (→ L2 halbe_exp) |
| 1942 | C. L. Siegel | Iteration of analytic functions, Ann. Math. 43 — `siegel1942` | Linearisierung für diophantische Drehzahlen | KAM-artige Abschätzung | bewiesen | G | fp.satz.neutral |
| 1944 | J. Hadamard; U. T. Bödewadt | Two works on iteration…, Bull. AMS 50; Zur Iteration reeller Funktionen, Math. Z. 49 — `hadamard1944`, `boedewadt1944` | reelle gebrochene Iteration: unendlich viele stetige Lösungen, Auswahlfrage | reelle Analysis | bewiesen | G | fg.kor.nichteindeutigkeit (nur eigen) |
| 1947 | E. Jabotinsky | C. R. 224 — `jabotinsky1947` | Komposition als Matrixprodukt; Anwendung auf $e^z$, $e^z-1$ | Bell-/Carleman-Matrizen | formal | G | — (→ L2 Carleman) |
| 1950 | R. Isaacs | Iterates of fractional order, Canad. J. Math. 2 — `isaacs1950` | Konstruktion gebrochener Iterierter reeller Funktionen | reell | bewiesen | G | — |
| 1955 | M. K. Fort | Embedding of homeomorphisms in flows, PAMS 6 — `fort1955` | wann ein Homöomorphismus in einen Fluss eingebettet werden kann | Topologie | bewiesen | G | — |
| 1955/58 | I. N. Baker | Math. Ann. 129; Zusammensetzungen ganzer Funktionen, Math. Z. 69 — `baker1955`, `baker1958` | keine ganze Lösung von $f(f(z))=F(z)$ für viele ganze $F$; formale Halbiterierte von $e^z-1$ haben Radius 0 | Wachstum ganzer Funktionen, Potenzreihen | bewiesen | Z | it.satz.baker (jetzt zitiert) |
| 1957 | S. Sternberg | Local contractions…, Amer. J. Math. 79 — `sternberg1957` | glatte ($C^\infty$) Linearisierung | Analysis | bewiesen | G | — |
| 1958 | G. Szekeres | Regular iteration of real and complex functions, Acta Math. 100 — `szekeres1958` | reguläre Iteration reeller Funktionen auch bei Multiplikator 1; durch Asymptotik eindeutig | Abel-Funktion als Grenzwert | bewiesen | G | it.satz.szekeres |
| 1958/60 | N. G. de Bruijn; W. J. Thron | Asymptotic Methods, Kap. 8; Sequences generated by iteration, Trans. AMS 96 — `debruijn1958`, `thron1960` | Asymptotik $x_n\sim(kan)^{-1/k}$ samt Log-Korrekturen nahe Multiplikator 1 | asymptotische Analysis | bewiesen | G | **—** (benutzt in fp.beob.tpid6_rate, fp.satz.p2_raten, kl.satz.eta_digamma) |
| 1960 | P. Erdős, E. Jabotinsky | On analytic iteration, J. Anal. Math. 8 — `erdos_jabotinsky1960` | Menge $S$ der Ordnungen mit analytischer Iterierter: ganz ℂ oder Mass 0 | formale Potenzreihen | bewiesen | R | it.satz.erdos_jabotinsky |
| 1962 | I. N. Baker | Permutable power series and regular iteration, J. Austral. MS 2 — `baker1962` | Gitteraussage: für $e^z-1$ konvergiert die Reihe nur für ganze Ordnungen | Vertauschbarkeit | bewiesen | Z | it.satz.baker, it.satz.erdos_jabotinsky |
| 1962–64 | G. Szekeres (mit K. W. Morris) | Fractional iteration of exponentially growing functions, J. Austral. MS 2; Tafeln ebd.; …entire and rational functions, 4 — `szekeres1961`, `szekeres1962_tafeln`, `szekeres1964` | Abel-Funktion an ∞ für schnell wachsende Funktionen; total monotone Abel-Funktion von $e^x-1$; Tafeln | Lévy-Methode verfeinert | bewiesen | Z (Hyperops-Wiki) | it.satz.szekeres_monoton |
| 1963 | E. Jabotinsky | Analytic iteration, Trans. AMS 108 — `jabotinsky1963` | Iterationsmatrizen, formale Iterationsgruppen | formal | bewiesen | G | it.satz.erdos_jabotinsky |
| 1964/67 | I. N. Baker | Fractional iteration near a fixpoint of multiplier 1; Non-embeddable functions…, Math. Z. 99 — `baker1964`, `baker1967` | bei Multiplikator 1 sind meromorphe Funktionen nur gebrochen linear in eine Gruppe einbettbar (laut Karlin–McGregor-Einleitung) | Wachstum | bewiesen | Z | — |
| 1965 | M. Lewin | Example of a function with non-analytic iterates, J. Austral. MS 5 — `lewin1965` | $e^z-1$ als Musterbeispiel nicht analytischer Iterierter | | bewiesen | Z | it.satz.baker |
| 1968 | S. Karlin, J. McGregor | Embedding iterates of analytic functions with two fixed points…, Trans. AMS 132 — `karlin_mcgregor1968` | **zwei Fixpunkte:** die regulären Iterationen beider Fixpunkte stimmen nur für gebrochen lineare $f$ überein | Koenigs an beiden Enden, Analytizität | bewiesen | Z (AMS-Einleitung, Forum tid 856) | it.satz.moebius_iteration, kl.beob.zwei_abel (jetzt zitiert); eigener Knoten fehlt |
| 1973–76 | J. Écalle | Théorie itérative, J. Math. Pures Appl. 54 (1975); C. R.-Noten 1971–76 — `ecalle1975` | holomorphe Invarianten parabolischer Keime; Gruppe der Iterationsordnungen; Resurgenz | Resurgenz, Borel-Summation | bewiesen | R (Bib.) | it.satz.ecalle_voronin |
| 1979 | M. C. Zdun | Continuous and differentiable iteration semigroups — `zdun1979` | reelle stetige/differenzierbare Iterationshalbgruppen | reell | bewiesen | G | — |
| 1980 | Rice, Schweizer, Sklar | When is f(f(z)) = az²+bz+c?, AMM 87 — `rice_schweizer_sklar1980` | funktionale Quadratwurzeln quadratischer Polynome | | bewiesen | G | — |
| 1981 | S. M. Voronin | Funct. Anal. Appl. 15 — `voronin1981` | analytische Klassifikation parabolischer Keime | Quasikonform | bewiesen | G | it.satz.ecalle_voronin |
| 1971–95 | Brjuno; Yoccoz; Herman 1985 | `brjuno1971`, `yoccoz1995`, `herman1985` | Brjuno-Bedingung hinreichend (Brjuno), für quadratische Polynome notwendig (Yoccoz); kritische Punkte auf Siegel-Rändern (Herman) | | bewiesen | G / Z | fp.satz.neutral, fp.satz.siegel_unbeschraenkt |
| 1989/2000 | P. Lavaurs; M. Shishikura | parabolische Implosion, gestörte Fatou-Koordinaten — `lavaurs1989`, `shishikura2000` | Grenzverhalten beim Aufbrechen parabolischer Punkte | | bewiesen | G | br.* (L3) |
| 1997 | S. C. Woon | Analytic continuation of operators, arXiv hep-th/9707206 — `woon1997` | „$s$-mal anwenden“ über die binomiale (Newton-)Reihe | Newton-Reihe | formal | A | kl.satz.newton_reihe (nur Forum) |
| 1998 | Aldrovandi, Freitas | Continuous iteration of dynamical maps, J. Math. Phys. 39 — `aldrovandi_freitas1998` | stetige Iteration über Bell-Matrizen | Matrizen | formal/numerisch | G | — (→ L2) |
| 2009/11 | Curtright, Zachos (Jin) | Evolution profiles…, J. Phys. A 42; Approximate solutions…, 44 — `curtright_zachos2009`, `curtright_jin_zachos2011` | Schröder-Methoden in der Physik, explizite Superfunktionen | Schröder | bewiesen/numerisch | Z | it.satz.elementare_superfunktionen |
| 2013 | G. A. Edgar | Fractional iteration of series and transseries, Trans. AMS 365 — `edgar2013` | reelle Iterationsfamilien für Transreihen (flach/moderat ja, tief nein) | Transreihen | bewiesen | R | — (→ L2 halbe_exp) |
| Lehrbücher | Kuczma 1968; Targonski 1981/1995; Kuczma–Choczewski–Ger 1990; Carleson–Gamelin 1993; Milnor 2006 | `kuczma1968`, `targonski1981`, `targonski1995`, `kuczma1990`, `carleson_gamelin1993`, `milnor2006` | Standarddarstellungen | | | G | vielfach |

## 2. Dieselbe Frage, verschiedene Zugänge

**F1 — Gibt es eine (eindeutige) analytische Iteration $f^t$ am Fixpunkt?** Hängt allein vom Multiplikator ab:
- hyperbolisch ($0<|\lambda|\ne1$): Schröder/Koenigs — holomorphe Linearisierung, Iteration eindeutig unter Holomorphie (it.satz.koenigs, it.satz.vertauschbar; Fatou 1924 für vertauschbare Abbildungen);
- superanziehend: Böttcher — Lücke;
- parabolisch ($\lambda=1$): drei Zugänge — (a) Leau/Fatou geometrisch (Blätter, Fatou-Koordinaten), (b) formal/algebraisch (Baker 1958/62, Erdős–Jabotinsky 1960, Jabotinsky 1963: formale Iterierte divergieren), (c) Klassifikation (Écalle über Resurgenz, Voronin quasikonform — dieselbe Antwort auf zwei Wegen); reell: Szekeres 1958 (Eindeutigkeit über Asymptotik);
- irrational indifferent: Cremer (nein, generisch), Siegel → Brjuno → Yoccoz (ja genau unter Brjuno, für quadratische Polynome).

**F2 — Gebrochene Iteration schnell wachsender Funktionen ($e^x$).** Lévy 1927/28 (über $e^x-1$, „reguläres Wachstum“), Kneser 1950 (über den komplexen Fixpunkt, L2), Szekeres 1962/64 (Abel-Funktion an ∞, total monoton), Jabotinsky 1947 / Walker 1991 (Matrizen, L2), Écalle/Edgar 2013 (Transreihen). Unterschied: Lévy/Szekeres liefern reelle, an ∞ reguläre Lösungen; Kneser die komplex-analytische. Ob sie übereinstimmen, ist eine Frage für L2/kl.

**F3 — Zwei Fixpunkte.** Satz: Karlin–McGregor 1968 (gleich nur bei Möbius). Numerisch: Asimov–Hickerson–Schroeppel um 1991, Forum 2007 („Bummer!“), Kouznetsov–Trappmann 2010. Eigen: Amplitudengesetz kl.verm.streifenbreite.

**F4 — Wie schnell konvergiert eine Iteration am parabolischen Punkt?** Leau/Fatou (Rate $n^{-1/k}$), Thron 1960, de Bruijn 1958 (mit Log-Korrekturen), Pomeau–Manneville 1980 (Engpass knapp neben dem Sattel-Knoten). Eigene Anwendungen: fp.beob.tpid6_rate, fp.satz.p2_raten, kl.satz.eta_digamma — alle drei sind Anwendungen dieser klassischen Rechnung.

## 3. Lücken (Priorität)

1. ~~hoch~~ erledigt 12.10. (it.satz.karlin_mcgregor) — Karlin–McGregor als eigener Literaturknoten (it), Kanten zu it.satz.moebius_iteration und kl.beob.zwei_abel.
2. ~~hoch~~ erledigt (it.satz.parabolische_rate) — „Rate am parabolischen Fixpunkt“ (Leau, Thron 1960, de Bruijn) als Literaturknoten; Kanten zu fp.beob.tpid6_rate, fp.satz.p2_raten, kl.satz.eta_digamma, fp.satz.parabolisch_rand.
3. ~~hoch~~ erledigt (it.satz.boettcher) — Böttcher-Satz (superanziehend) fehlt ganz.
4. erledigt (it.konstr.levy) — Lévys Weg (1927/28) als Knoten, abgestimmt mit L2 (halbe Exponentialfunktion).
5. teils erledigt (it.satz.cremer; Quellen an it.satz.leau_fatou, it.satz.vertauschbar ergänzt) — Cremer-Punkte als eigener Knoten bzw. Teil von fp.satz.neutral; Leau 1897 und Fatou 1919 als Quellen an it.satz.leau_fatou; Fatou 1924 an it.satz.vertauschbar.
6. Knoten it.satz.baker_einbettung (nur Wiedergabe geprüft) — Baker 1964/1967 (Nicht-Einbettbarkeit bei Multiplikator 1) im Original prüfen.
7. **gering** — Sternberg, Fort, Zdun (glatte bzw. reelle Einbettung), Babbage, Rice–Schweizer–Sklar, Curtright–Zachos, Aldrovandi–Freitas.
8. Nicht zugänglich (Volltext): Lévy 1928, Baker 1958/62, Szekeres 1958 — Aussagen nur über Rezensionen und Zitate.

## 4. Stand und Restlücken (Abschluss L1, 12.10.2026)

- **Stand:** 12 → 19 Aussagen (it.*). Neu als Literaturknoten (`daten/aussagen/lit_it_1.yaml`): it.satz.karlin_mcgregor, it.satz.parabolische_rate, it.satz.boettcher,
  it.konstr.levy, it.satz.cremer, it.satz.baker_einbettung, it.satz.sternberg_reell. Quellen nachgetragen bei Koenigs/Leau–Fatou/Baker/Erdős–Jabotinsky/vertauschbar/Möbius/
  elementare Superfunktionen. Neuer Strang it.indifferent.
- **Restlücken:** Lévy 1928, Baker 1958/62/64/67, Szekeres 1958 und Écalle nur über Rezensionen bzw. Zitate (kein Volltext gelesen); Fatou 1924 als eigener
  Knoten (vertauschbare Abbildungen); Jabotinsky 1947 / Aldrovandi–Freitas (Matrizen) und Edgar 2013 (Transreihen) liegen bei L2; Rice–Schweizer–Sklar 1980 und
  Babbage nur in der Karte. Fort 1955 und Zdun 1979: genauer Satz nicht gelesen (nur als Literatur bei fg.satz.einbettung_fluss).

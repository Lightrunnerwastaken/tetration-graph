# Literaturkarte kl — Lösungsklassen nach Basis (L1, Stand 12.10.2026)

Prüfvermerk wie in KARTE_it.md (O Original · R Rezension · A Abstract · Z Zitat/Wikipedia/Forum · G Gedächtnis).

Bestand vorher: 36 Aussagen (kl.*) in 04_ und 10_ (Bruchsummen/Basiswechsel); Literatur vor allem Kouznetsov–Trappmann
2010/2012, Müller–Schleicher, Forumsfäden. Viel Eigenes (Bruchsummen bei η, Basiswechsel an Fixpunkten), dazu mehrere „neu“-Einstufungen.

## 1. Welche Klasse für welche Basis — wer hat was gezeigt

| Basisbereich | Dynamik von $b^x$ | Klassische Lösung(en) | Wer | Beweisstand | Bei uns |
|---|---|---|---|---|---|
| $0<b<e^{-e}$ | anziehender 2-Zyklus | Turm pendelt (Barrow 1936, Knoebel 1981); Tetration als Funktion kaum behandelt | Barrow; Baker–Rippon 1985; Galidakis 2005 | Turm bewiesen; Tetration offen | fp.satz.p2_form, fp.satz.rechteckwelle (eigen) |
| $e^{-e}\le b<1$ | anziehender Fixpunkt, $\lambda<0$ | reguläre Iteration mit $\lambda^t$ komplex — keine reelle reguläre Tetration | (keine Arbeit gefunden; evtl. Paulsen 2019) | — | **Lücke** |
| $1<b<\eta$ | zwei reelle Fixpunkte $a<A$ | reguläre Tetration am unteren Fixpunkt (Koenigs); vier reguläre Superexponentiale (je zwei an $a$ und $A$) | Koenigs 1884; Kouznetsov–Trappmann 2010; Karlin–McGregor 1968 (an $a$ und $A$ verschieden) | bewiesen | kl.satz.regulaer, kl.satz.vier_sqrt2, kl.beob.zwei_abel, kl.satz.obere_superexp |
| $b=\eta$ | parabolischer Fixpunkt $e$ | Fatou-Koordinaten der zwei Blätter: zwei reguläre Superexponentiale; „cheta“ von oben | Leau/Fatou; Écalle; Trappmann–Kouznetsov 2012; Fox 2009 (Forum) | bewiesen | kl.satz.parabolisch, kl.satz.cheta, kl.kor.eta_baker |
| Übergang an η | | ob die Tetration in $b$ analytisch über η geht | Kouznetsov (Mizugadro, Taylor-Tafel), mike3 2011 (Forum, eher Verzweigung) | offen | kl.beob.stetig_eta, kl.frage.stetig_eta, kl.frage.eta_analytisch |
| $b>\eta$ | nur komplexe Fixpunkte | Kneser; Eindeutigkeit Paulsen–Cowgill | Kneser 1950; Kouznetsov 2009; Paulsen–Cowgill 2017; Trappmann–Kouznetsov 2011 | bewiesen | kn.* (L2) |
| komplexe $b$ | Shell–Thron innen: anziehend | reguläre Tetration; holomorph in $b$ | Paulsen 2019 | bewiesen | kn.satz.paulsen_komplex (L2), fp.satz.shell_thron |

## 2. Werke zu Auswahl- und Eindeutigkeitsprinzipien

| Jahr | Autor(en) | Werk, Ort — Schlüssel | Hauptergebnis | Methode | Beweisstand | Prüf. | Bei uns |
|---|---|---|---|---|---|---|---|
| 1914 | F. Carlson | Sur une classe de séries de Taylor, Diss. Uppsala (Schlüssel folgt: `carlson1914`) | ganze Funktion vom Exponentialtyp, Typ $<\pi$ auf der imaginären Achse, Nullstellen in ℕ ⇒ $\equiv0$ | Phragmén–Lindelöf | bewiesen | Z (Wikipedia) | kl.satz.carlson (nur Forum tid 496) |
| 1924 | N. E. Nörlund | Vorlesungen über Differenzenrechnung (Bib?) | Konvergenz von Newton-Reihen in Halbebenen | | bewiesen | G | kl.satz.newton_reihe |
| 1927/28 | P. Lévy | `levy1927`, `levy1928` | „reguläres Wachstum“ als Auswahlprinzip für gebrochene Iteration | | bewiesen | R | — |
| 1961/62 | G. Szekeres | `szekeres1961` | total monotone Abel-Funktion von $e^x-1$ — reelle Eindeutigkeitsbedingung | | bewiesen | Z | it.satz.szekeres_monoton |
| 1968 | S. Karlin, J. McGregor | `karlin_mcgregor1968` | reguläre Iterationen an zwei Fixpunkten gleich ⇔ Möbius | | bewiesen | Z | kl.beob.zwei_abel (jetzt zitiert) |
| 1991 | P. Walker | `walker1991` | $C^\infty$-Lösungen über Abel-Funktion von $e^x-1$ / Matrizen | | bewiesen | Z | kn.konstr.walker (L2) |
| ~1991 | Asimov, Hickerson, Schroeppel | unveröffentlicht, über Forum tid 69 | reguläre Lösungen an 2 und 4 zur Basis $\sqrt2$ verschieden ($\sim10^{-24}$) | numerisch | numerisch | Z | kl.beob.zwei_abel |
| 1997 | S. C. Woon | `woon1997` | Newton-/Binomialreihe für „$s$-mal anwenden“ | | formal | A | kl.satz.newton_reihe (nur Forum) |
| 2005/11 | M. Müller, D. Schleicher | `mueller_schleicher2005`, `mueller_schleicher2011` | Bruchsummen axiomatisch eindeutig | | bewiesen | Z | kl.satz.bruchsumme, kl.satz.bruchprodukt, kl.satz.eta_digamma |
| 2006 | M. H. Hooshmand | `hooshmand2006` | Eindeutigkeit unter Konvexität | | bewiesen | Z | kn.satz.hooshmand (L2) |
| 2010 | D. Kouznetsov, H. Trappmann | Portrait of the four regular super-exponentials to base √2, Math. Comp. 79 — `kouznetsov_trappmann2010` | Begriff reguläre Superfunktion; Bedingungen für ganze Superfunktionen; zwei Eindeutigkeitskriterien; vier reelle reguläre Superexponentiale | Koenigs, Cauchy-Integral | bewiesen | A | kl.satz.regulaer, kl.satz.vier_sqrt2 |
| 2010 | Forum (tommy1729; R. Israel über sci.math) | tid 476 — `forum_helms2010` | $e-\eta\uparrow\uparrow k\approx2e/k$, Summe divergiert logarithmisch | heuristisch | heuristisch | Z | kl.satz.eta_digamma (jetzt korrigiert) |
| 2012 | H. Trappmann, D. Kouznetsov | Computation of the two regular super-exponentials to base $e^{1/e}$, Math. Comp. 81 — `kouznetsov_trappmann2012` | zwei reguläre Superexponentiale bei η | Fatou-Koordinaten, Reihen | bewiesen | A | kl.satz.parabolisch, kl.satz.eta_koordinate |
| 2011/17 | Paulsen–Cowgill; Trappmann–Kouznetsov | `paulsen_cowgill2017`, `trappmann_kouznetsov2011` | Eindeutigkeit über ganzes Streifenverhalten | | bewiesen | Z | kn.* (L2) |
| 2015/21 | J. D. Nixon | `nixon2015`, `nixon2021a` | beschränkte analytische Hyperoperatoren; Frage nach $\mathrm{Re}\,T\ge1$ in der rechten Halbebene | | teilweise | G | kl.satz.folgen_vollmonoton |
| 2017 | V. Reshetnikov | MO 259278 — `reshetnikov2017` | q-Binomial-Reihe der Tetration unter η; Vermutung totaler Monotonie | | numerisch | O (laut Bestand) | kl.beob.q_binomial, kl.verm.vollmonoton |

## 3. Dieselbe Frage, verschiedene Zugänge

**F1 — Welche Lösung ist „die richtige“ unter η?** (a) Holomorphie am unteren Fixpunkt (Koenigs, Kouznetsov–Trappmann: reguläre Tetration); (b) Beschränktheit in der rechten Halbebene + Carlson (Forum tid 496; kl.satz.carlson); (c) Newton-Interpolation der ganzen Werte (Woon, Forum; kl.satz.newton_reihe); (d) totale Monotonie / Vorzeichen der Koenigs-Koeffizienten (Reshetnikov, Nixon; eigene Bedingung kl.satz.folgen_vollmonoton); (e) Bruchsummen-Axiome (Müller–Schleicher, kl.satz.bruchprodukt). Bei uns sind (a)–(e) als gleichwertig für die reguläre Tetration gezeigt oder belegt — das ist die eigene Synthese.

**F2 — Was unterscheidet die beiden Fixpunkte?** Satz (Karlin–McGregor), Zahl (Asimov–Hickerson–Schroeppel/Forum, Kouznetsov–Trappmann), eigenes Amplitudengesetz (kl.verm.streifenbreite, kl.verm.faktor).

**F3 — Was passiert bei η?** Écalle/Fatou (parabolisch), Trappmann–Kouznetsov 2012 (zwei Lösungen), Forum (cheta, Fox 2009); Stetigkeit in $b$: Kouznetsov (Taylor-Tafel) gegen mike3 (Verzweigung) — eigene Messung kl.beob.stetig_eta ($C^1$ numerisch).

## 4. Befunde zur Neuheit (12.10.2026)

- **kl.satz.eta_digamma**: Divergenz $\sim2e\ln n$ und $e-\eta\uparrow\uparrow k\approx2e/k$ (sogar mit $\ln k/3$-Korrektur) standen schon 2010 im Forumsthread tid 476 (tommy1729, R. Israel). Nur Digamma-Form, Rest $R_T$ und Beweis sind eigen. → Recherche-Text korrigiert.
- **kl.beob.zwei_abel**: Nicht-Konstanz folgt aus Karlin–McGregor 1968. → Prüffrage f1 und Quelle ergänzt.
- **kl.verm.streifenbreite / kl.verm.faktor**: Volltext Kouznetsov–Trappmann 2010 nicht gelesen; Abstract sagt nichts zum Fourier-Abfall. Bleibt „vermutlich neu“, Lücke 2.
- **kl.satz.bruchprodukt, kl.satz.moebius_qdigamma, kl.satz.basiswechsel_gleiches_L, kl.beob.eta_konstanten**: nichts Vorwegnehmendes gefunden (nur Abstract-Ebene).

## 5. Lücken (Priorität)

1. erledigt 12.10. (Quellen an kl.satz.carlson, kl.satz.newton_reihe) — Carlson 1914 (und Hardy 1920, Rubel 1956) als Literaturquelle an kl.satz.carlson; Woon 1997 an kl.satz.newton_reihe.
2. **hoch** — Volltext Kouznetsov–Trappmann 2010 auf Aussagen zum Fourier-Abfall der Differenz prüfen (kl.verm.streifenbreite).
3. **mittel** — Basen $e^{-e}\le b<1$ und $0<b<e^{-e}$: Literatur zur Tetration (nicht nur zum Turm) suchen; Paulsen 2019 prüfen (gemeinsam mit L2).
4. **mittel** — Kouznetsov–Trappmann 2010: die zwei Eindeutigkeitskriterien als eigene Knoten (mit L2 abstimmen, falls schon kn.satz.trappmann_kouznetsov).
5. **mittel** — Lévys Auswahl durch „reguläres Wachstum“ als Zugang (F1) — mit it-Lücke 4.
6. **gering** — Nörlund für die Konvergenz der Newton-Reihe.

## 6. Stand und Restlücken (Abschluss L1, 12.10.2026)

- **Stand:** keine neuen kl-Knoten, aber Neuheit und Quellen korrigiert: kl.satz.eta_digamma (Forum 2010 hatte 2e/k), kl.beob.zwei_abel (Karlin–McGregor),
  kl.kor.eta_baker (Baker 1958/62), kl.satz.carlson (Carlson 1914, Hardy 1920), kl.satz.newton_reihe (Woon 1997), kl.verm.streifenbreite (Forum tid 69 ohne Formel).
- **Restlücken:** Volltext Kouznetsov–Trappmann 2010 (Fourier-Abfall der Differenz; frei nur als PDF — Leitung/Autor); Basen $e^{-e}\le b<1$: keine Arbeit zur
  Tetration (nur zum Turm) gefunden — Forum tid 352 (2009) stellt die Frage für $b<e^{-e}$, ungelöst; Paulsen 2019 auf Basen unter 1 prüfen (L2);
  Kouznetsov–Trappmanns zwei Eindeutigkeitskriterien als eigene Knoten; Nörlund (Newton-Reihen).

**Nachtrag 12.10. 18:15:** Kouznetsov–Trappmann 2010 (freies Manuskript) gelesen: §5, (5.9)–(5.10) — erste Mode der Differenz mit
$\xi=0.7252\cdot10^{-24}$, gültig mindestens für $\mathrm{Im}\,z<T_2/2$. Der Streifen ist dort gesehen, das Abfallgesetz der höheren Moden
(kl.verm.streifenbreite, kl.verm.faktor) nicht — beide bleiben „vermutlich neu“. Lévy 1928 nur hinter Bezahlschranke (Springer) — Restlücke.

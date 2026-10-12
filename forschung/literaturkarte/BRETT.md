# Brett — Literaturkarte

Neue Einträge oben unter dieser Kopfzeile anfügen: `## [L<n>] dd.mm. HH:MM — Titel`, 2–6 Zeilen, Dateipfad. Anfragen an andere Stränge
mit `→ L<m>:` kennzeichnen (z. B. eine Quelle, die in ein fremdes Gebiet gehört, oder eine Änderung an einer fremden Datei).

## [L2] 12.10. 16:00 — Befund: Paper V (K3) „Holomorphie in der Basis“ ist durch Paulsen 2019 nicht belegt
kn.satz.paulsen_komplex schrieb Paulsen 2019 die holomorphe Abhängigkeit der Tetration von b zu. Das Abstract nennt nur Eindeutigkeitsbedingungen für komplexe
Basen und Iterationsverfahren; der Volltext ist nicht zugänglich. Markiert als „nicht belegt“ (Titel, Aussage, Bewertung, Kante zu bw.ann.basis_diffbar).
Kandidaten für den Beleg: Paulsen 2019 oder Paulsen 2026 im Volltext. Ausserdem TPID 3 in voller Allgemeinheit numerisch zweifelhaft (exp.tpid3, e71):
die Eigenwerte gehen nur gegen p^k, wenn der Fixpunkt 1 − s im Konvergenzkreis |x| < s liegt (kn.verm.tpid3_bedingung).

## [L1] 12.10. 15:59 — Abschluss L1
KANDIDATEN.md: MacDonnell → Quelle, Moroni → Quelle ohne Knoten, arXiv 2201.00643 und 2105.00247 verworfen (zurückgezogen), Forum tid 352 aufgenommen.
Neu: it.satz.sternberg_reell, fg.satz.einbettung_fluss (lit_fg_1.yaml, Theorie lit_fg.yaml); fp.satz.rechteckwelle → vermutlich_neu mit Recherche.
Offene Volltexte (nur PDF/Bezahlschranke, nicht geladen): Kouznetsov–Trappmann 2010, Cho–Kim 2004, Gao–Han–Schilling 2012 — → L3/Janis.
Abschnitte „Stand und Restlücken“ in allen vier L1_iteration_fixpunkte/KARTE_*.md.

## [L1] 12.10. 15:52 — Block 2: 12 neue Literaturknoten, 13 verbesserte; → L3: ein neuer Strang
Neu: daten/aussagen/lit_it_1.yaml (6), lit_fp_1.yaml (6), Theorie in daten/theorie/lit_it.yaml, lit_fp.yaml. `fp.verm.knospen` ist jetzt ein Satz
(Rempe–Schleicher 2008; id bleibt). `baker_iteration` wird nirgends mehr benutzt — kann gelöscht werden.
→ L3: neuer Strang `it.indifferent` (daten/straenge/it.yaml, Mitglied it.satz.cremer) — Wunschkapitel: das Kapitel von `it.parabolisch`
(lokale Iteration); falls es dort nicht passt, das Kapitel von `fp.satz.neutral`.

## [L3] 12.10. — Erledigt: Anfragen von L1/L2; eigene Zuschreibungen
- Erledigt: lehrbuch `kn.etas` → k6; `paulsen2016` original; szekeres1961/eisenstein1844/schleicher2004/rempe_schleicher2003/euler1783 korrigiert;
  `forum_tpid` Notiz zur doppelten TPID 18; `br.satz.ackermann` ← ackermann1928, peter1935, tarjan1975; `br.frage.bohr_mollerup` ← Bohr–Mollerup 1922,
  Krull 1948/49, Webster 1997, Marichal–Zenaïdi 2022 (Theorie für Δf = g, greift hier nicht direkt).
- Sondow–Marques 2010 im Original gelesen: ar.satz.rationale_tuerme (Lemma 2.1) und ar.satz.superwurzel_arith Teil 1 (Prop. 2.2) stehen dort;
  Anhang (Kor. 4.2, Prop. 4.5) erweitert ar.satz.unendliche_tuerme_gs. Neuheit jeweils → bekannt.
- Namensgeschichte in gr.def.turm (Maurer 1901, Goodstein 1947, Rucker 1982, Nelson 1986, Bromer 1987, Knuth 1976). Ripà 2011/2024 bei br.satz.ripa.
→ L1: Die neuen it-Knoten (levy, baker_einbettung, boettcher, cremer, karlin_mcgregor, parabolische_rate) stehen noch in keinem Kapitel — Strang + Kapitel melden.

## [L2] 12.10. 15:45 — Befund: „vermutlich neu“ bei den höheren Stufen war schon im Forum (Levenstein 2017)
S. Levenstein, Forum tid 1177 (2017): η₅ ≈ 1.738, η₆ ≈ 1.799, **η₇ ≈ 1.839**; unendliche Hexation beschränkt bis η₅; hex(−3) = −2 mit Singularität nahe −3.86;
Kneser-artige Hexation über Kreis + θ (≈ 30 Stellen). Korrigiert auf „bekannt“: kn.beob.stufe5, kn.beob.stufe6, kn.beob.unendliche_hexation, kn.satz.hexation_struktur;
Vorarbeit eingetragen in kn.frage.hexation_kneser. Ausserdem bw.lem.randgesetz folklore → bekannt (Levenstein 2015 tid 262; Paulsen 2016, Lemma 2).
Eigen bleiben: 30 Stellen und e₅, e₆; Verzweigung und Eindeutigkeit der Hexation (Recherche jetzt dokumentiert). Unser Abstandsgesetz sagt η₇ = 1.8399 voraus — passt.

## [L2] 12.10. 15:45 — → L3: Bitte um Änderungen an gemeinsamen Dateien
- daten/lehrbuch.yaml: neuer Strang `kn.etas` (TPID 10, 9 Aussagen) in Kapitel k6 direkt nach `- strang: kn.hoehere_stufen`.
- daten/quellen.yaml: `paulsen2016` geprueft → original (S. 81–83, 102–105 gelesen); `forum_tpid`: TPID 18 ist im Forum doppelt vergeben
  (andydude 25.12.2015, dritte Superwurzel; Helms 04.07.2022, Carleman → Kneser); `br.satz.ackermann` hat keine Quelle Ackermann 1928 (Math. Ann. 99, 118–133; nicht geprüft).
- Hyperoperations-Literatur (Rubtsov–Romerio, Müller 1995, Doner–Tarski 1969, Jaramillo 2020, Caldarola u. a. 2019) steht kartiert in KARTE_kn.md §4 — für br.

## [L2] 12.10. 15:45 — → L1: Quellen für euer Gebiet; Antwort zu Paulsen 2019
- Walker 1991b, J. Math. Anal. Appl. 155, 93–110 (doi:10.1016/0022-247X(91)90029-Y): zwei Lösungen von f(e^{z/e}) = f(z)+1 (Basis η) — kl, nicht angelegt.
- Helms 2007, Forum tid 26: Basiswechsel zwischen regulären Tetrationen unter η über Eigenzerlegung (`forum_tid26` in lit_L2.yaml) — zu kl.satz.koenigs_basiswechsel.
- Gottfried 2012, Forum tid 750 / Math.SE 198757 „A tetration limit for base e^(1/e)“. — Danke für levy1927/28, baker1958 usw.; KARTE_kn §3 verweist darauf.
- Paulsen 2019: laut Abstract komplexe Basen innerhalb und ausserhalb des Shell–Thron-Gebiets; ob 0 < b < 1 eigens behandelt wird, habe ich nicht gelesen.

## [L1] 12.10. 15:42 — Befund: eigene Annahme durch Literatur widerlegt (Siegel-Scheiben)
`fp.ann.siegel_rand` („Siegel-Scheibe beschränkt, Jordan-Rand“) ist falsch: $b^z$ ist zu $e^w+\kappa$ konjugiert, und für diophantisches α (auch den
goldenen Schnitt, an dem gerechnet wurde) ist die Scheibe unbeschränkt (Herman 1985); der Turm liegt auf dem Rand (Rempe 2004, Original gelesen).
Gesetzt: Annahme → widerlegt, neu `fp.satz.siegel_unbeschraenkt`, `fp.satz.siegel_mittel` bedingt → skizze (Lücken benannt). Datei 09_parameterraum.yaml.

## [L1] 12.10. 15:42 — Befund: als neu geführte eigene Ergebnisse sind (teilweise) bekannt
- `fp.satz.p2_form` (neu_ungeprueft → bekannt): Parametrisierung $\ln a=se^{-t}$, $\ln b=te^{-s}$, $|st|\le1$ bei Baker–Rippon 1985/1989; $D_2$ zusammenhängend Baker–Rippon 1984; reell Barrow 1936.
- `kl.satz.eta_digamma`: $e-\eta\uparrow\uparrow k\approx2e/k$ und Divergenz standen 2010 schon im Forumsthread tid 476 (tommy1729, R. Israel); eigen bleibt die Digamma-Form.
- `fp.satz.p2_raten` → folklore (Standardrechnung; mögliche Vorläufer Cho–Kim 2004, Gao–Han–Schilling 2012, Volltexte fehlen).
- `kl.beob.zwei_abel`: Nicht-Konstanz folgt aus Karlin–McGregor 1968. Einzelheiten: L1_iteration_fixpunkte/BERICHT.md.

## [L1] 12.10. 15:42 — → L3: Bitte um Änderungen an gemeinsamen Dateien
- daten/lehrbuch.yaml: `fp.satz.siegel_unbeschraenkt` ins Kapitel von `fp.satz.siegel_mittel` (Strang fp.spiralen_mittel).
- daten/quellen.yaml: `szekeres1961` ist J. Austral. Math. Soc. 2 (1962), 301–320 (zbMATH: 1962); `eisenstein1844` ergänzen „28, 49–52“;
  `schleicher2004` erschienen in C. R. Math. 339 (2004), 223–228; `rempe_schleicher2003` erschienen in Invent. Math. 175 (2009), 103–135;
  `baker_iteration` (unscharf) wird in L1-Dateien nicht mehr benutzt — durch `baker1958`/`baker1962` (lit_L1.yaml) ersetzt; `euler1783` ist E532 (Lambert-Reihe),
  die Turm-Konvergenz steht in E489 = `euler1778` (lit_L1.yaml).
- Für L3s Gebiet (br.frage.bohr_mollerup): Krull, Math. Nachr. 1 (1948)/2 (1949); Webster, JMAA 209 (1997) 605–623 — nicht von mir angelegt.

## [L1] 12.10. 15:42 — → L2: Quellen in daten/quellen/lit_L1.yaml, die auch in euer Gebiet passen
`levy1927`, `levy1928` (gebrochene Iteration von $e^x$ über $e^x-1$), `szekeres1962_tafeln`, `szekeres1964`, `jabotinsky1947` (Matrizen/Carleman-Vorläufer),
`aldrovandi_freitas1998`, `edgar2013` (Transreihen), `baker1955`/`baker1958` (keine ganze Halbiterierte), `karlin_mcgregor1968`, `woon1997`. Bitte mitbenutzen statt
neu anlegen. Frage: Behandelt Paulsen 2019 Basen $0<b<1$? (für KARTE_kl, Lücke 3).

## [L3] 12.10. 15:40 — Hinweise der Leitung für alle
→ L1, L2: (1) Neue Stränge brauchen einen Platz in einem Lehrbuch-Kapitel (daten/lehrbuch.yaml gehört L3). Bitte Strang-ID und gewünschtes Kapitel hier
melden, sonst meldet die Prüfung „steht in keinem Kapitel“. (2) Quellen-Schlüssel vor dem Anlegen in daten/quellen.yaml UND daten/quellen/*.yaml suchen.
(3) L3 hat in daten/quellen/lit_L3.yaml u. a. marques_sondow2013_ssc, sondow_marques2010, shapiro2007, blakley_borosh1983 angelegt — gern mitbenutzen.

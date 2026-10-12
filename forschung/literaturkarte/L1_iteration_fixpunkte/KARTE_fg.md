# Literaturkarte fg — Funktionalgleichungen (L1, Stand 12.10.2026)

Prüfvermerk wie in KARTE_it.md (O Original · R Rezension · A Abstract · Z Zitat/Wikipedia/Forum · G Gedächtnis; „Bib?“ = Bibliografie nicht nachgeprüft).

Bestand vorher: 12 Aussagen (fg.*): 4 Definitionen, 8 Sätze; als Literatur fast nur Kuczma–Choczewski–Ger 1990 (ohne Satznummern) und Koenigs.
Die Gleichungen selbst (Abel, Schröder, Böttcher, Julia) sind klassisch; eigene Arbeit steckt in der Anwendung auf $F(x+1)=b^{F(x)}$.

## 1. Werke

| Jahr | Autor(en) | Werk, Ort — Schlüssel | Hauptergebnis | Methode | Beweisstand | Prüf. | Bei uns |
|---|---|---|---|---|---|---|---|
| 1815 | C. Babbage | Phil. Trans. 105 — `babbage1815` | Involutionen $f(f(x))=x$, funktionale Wurzeln; Konjugation erzeugt Lösungen | Konjugation | Beispiele | Z | — |
| ~1824/1881 | N. H. Abel | Œuvres II — `abel1881` | Abel-Gleichung $\alpha\circ f=\alpha+1$ | | formal | G | fg.def.abel (nur Kuczma) |
| 1871 | E. Schröder | Math. Ann. 3 — `schroeder1871` | Schröder-Gleichung $\sigma\circ f=\lambda\sigma$ | | Beispiele | R, Z | fg.def.schroeder (Koenigs, Kuczma) |
| 1884 | G. Koenigs | Ann. ÉNS — `koenigs1884` | Existenz/Eindeutigkeit der holomorphen Schröder-Lösung | | bewiesen | R | fg.satz.schroeder_abel |
| 1904 | L. E. Böttcher | `boettcher1904` | Böttcher-Gleichung $\phi\circ f=\phi^k$ | | bewiesen | G | — |
| 1918/22 | G. Julia | `julia1918` | Julia-Gleichung $\lambda(f(x))=f'(x)\lambda(x)$ für die Erzeugende eines Iterationsflusses | | bewiesen | R (Bib.) | it.lem.julia (nur Kuczma) |
| 1948/49 | W. Krull | Bemerkungen zur Differenzengleichung $g(x+1)-g(x)=\varphi(x)$ I/II, Math. Nachr. 1/2 — (Schlüssel bei L3) | Eindeutigkeit konvexer/monotoner Lösungen von Differenzengleichungen (Bohr–Mollerup-Typ) | Konvexität | bewiesen | R (Teil II) | → L3 br.frage.bohr_mollerup |
| 1955 | M. K. Fort | PAMS 6 — `fort1955` | Einbettung von Homöomorphismen in Flüsse (reelle Iterationsgruppen) | Topologie | bewiesen | G | fg.satz.allgemeine_loesung (nur indirekt) |
| 1966 | J. Aczél | Lectures on Functional Equations and Their Applications, Academic Press (Bib?) | Standardwerk allgemeiner Funktionalgleichungen | | | G | — |
| 1968 | M. Kuczma | Functional Equations in a Single Variable — `kuczma1968` | allgemeine Lösung der Abel-/Schröder-Gleichung = spezielle Lösung ∘ (id + periodisch); Eindeutigkeit unter Monotonie/Konvexität/Regularität am Fixpunkt | reelle Analysis | bewiesen | G | fg.satz.allgemeine_loesung, fg.kor.nichteindeutigkeit |
| 1979 | M. C. Zdun | `zdun1979` | stetige und differenzierbare Iterationshalbgruppen | | bewiesen | G | — |
| 1981/95 | G. Targonski | Topics in Iteration Theory; Progress… since 1981 — `targonski1981`, `targonski1995` | Übersichten | | | G | — |
| 1990 | Kuczma, Choczewski, Ger | Iterative Functional Equations, CUP — `kuczma1990` | Standardmonografie (Abel, Schröder, Böttcher, Julia, Iterationsgruppen) | | bewiesen | G | fg.def.*, fg.satz.* |
| 1997 | R. Webster | Log-convex solutions to $f(x+1)=g(x)f(x)$, JMAA 209 — (L3) | Γ-artige Funktionen, Eindeutigkeit über Log-Konvexität | Konvexität | bewiesen | R | → L3 |
| 2001 | K. Baron, W. Jarczyk | Recent results on functional equations in a single variable, Aequationes 61 — `baron_jarczyk2001` | Übersicht seit 1990, offene Probleme, 281 Literaturangaben | | | R | — |
| 2003 | G. Belitskii, V. Tkachenko | One-dimensional Functional Equations, Birkhäuser — `belitskii_tkachenko2003` | Lösbarkeit linearer und kohomologischer Gleichungen $\varphi(f(x))-\varphi(x)=g(x)$ | | bewiesen | G | fg.satz.kohomologisch (nur Kuczma, eigen) |
| 2005/11 | M. Müller, D. Schleicher | Bruchsummen, J. Comput. Appl. Math. 178; AMM 118 — `mueller_schleicher2005`, `mueller_schleicher2011` | Differenzengleichung $S(x+1)-S(x)=f(x)$ mit Axiomen eindeutig | Asymptotik | bewiesen | Z | kl.satz.bruchsumme (L1/kl) |
| 2006 | M. H. Hooshmand | Ultra power and ultra exponential functions, ITSF 17 — `hooshmand2006` | Eindeutigkeit der Tetration unter Konvexitätsbedingung | Konvexität | bewiesen | Z | kn.satz.hooshmand (L2) |
| 2010 | D. Kouznetsov, H. Trappmann | Math. Comp. 79 — `kouznetsov_trappmann2010` | Begriff Superfunktion, ganze Superfunktionen, Eindeutigkeitskriterien | | bewiesen | A | fg.kor.keine_ganze_tetration, kl.* |
| 2014 | M. C. Zdun, P. Solarz | Aequationes 87 — `zdun_solarz2014` | Übersicht reelle Iterationsgruppen/-halbgruppen | | | G | — |

## 2. Dieselbe Frage, verschiedene Zugänge — „Wie viele Lösungen hat $F(x+1)=b^{F(x)}$, und was macht eine eindeutig?“

| Zugang | Wer | Auswahlbedingung | Ergebnis | Bei uns |
|---|---|---|---|---|
| reell-topologisch | Kuczma 1968, Fort 1955, Zdun 1979 | keine | alle stetigen monotonen Lösungen = $F\circ(\mathrm{id}+\theta)$, $\theta$ 1-periodisch: unendlich viele | fg.satz.allgemeine_loesung, fg.kor.nichteindeutigkeit |
| Konvexität | Krull 1948, Webster 1997, Hooshmand 2006 | (log-)konvex o. ä. | Eindeutigkeit wie bei Bohr–Mollerup | kn.satz.hooshmand (L2), br.frage.bohr_mollerup (L3) |
| Holomorphie am Fixpunkt | Koenigs 1884, Karlin–McGregor 1968 | $F$ regulär am Fixpunkt | reguläre Tetration (unter η), bei zwei Fixpunkten verschieden | kl.satz.regulaer, kl.beob.zwei_abel |
| komplexes Wachstum | Kneser 1950; Paulsen–Cowgill 2017; Trappmann–Kouznetsov 2011 | holomorph im Streifen, Grenzwerte bei $\pm i\infty$ | Kneser-Lösung eindeutig | kn.* (L2) |
| totale Monotonie | Szekeres 1961 ($e^x-1$); Reshetnikov 2017, Nixon (Tetration) | Abel-Funktion total monoton / Koenigs-Koeffizienten ≤ 0 | eindeutig bzw. vermutet | it.satz.szekeres_monoton, kl.verm.vollmonoton |
| Bruchsummen-Axiome | Müller–Schleicher 2011 | Asymptotik der Summanden | eindeutige Bruchsumme | kl.satz.bruchsumme |

## 3. Lücken (Priorität)

1. **mittel** — Abel 1881 und Schröder 1871 als Quellen an fg.def.abel / fg.def.schroeder (Ursprung der Gleichungen); Babbage als Geschichte.
2. **mittel** — Böttcher-Gleichung als Definition + Satz (zusammen mit it-Lücke 3).
3. **mittel** — Satznummern bei Kuczma 1968 / Kuczma–Choczewski–Ger für fg.satz.allgemeine_loesung, fg.satz.kohomologisch: braucht das Buch (Leitung fragen).
4. **mittel** — Belitskii–Tkachenko für fg.satz.kohomologisch; Julia-Gleichung mit Julia/Écalle-Quelle an it.lem.julia.
5. **gering** — Übersichten Baron–Jarczyk 2001, Zdun–Solarz 2014, Targonski als Einstieg ins Gebiet (Theorie-Block „einordnung“).
6. An L3: Krull 1948/49 und Webster 1997 gehören zu br.frage.bohr_mollerup.

## 4. Stand und Restlücken (Abschluss L1, 12.10.2026)

- **Stand:** 12 → 14 Aussagen (fg.*). Neu (Block 3): fg.def.boettcher; Fort 1955 jetzt sekundär belegt (Intervall-Homöomorphismus einbettbar ⇔ ordnungserhaltend, nach Jones 1972). Neu: fg.satz.einbettung_fluss (`daten/aussagen/lit_fg_1.yaml`, Beweis vollständig, Literatur Kuczma/Fort/Zdun).
  Ursprungsquellen Abel 1881 und Schröder 1871 an fg.def.abel / fg.def.schroeder.
- **Restlücken:** Satznummern bei Kuczma 1968 und Kuczma–Choczewski–Ger 1990 (Buch nötig); Belitskii–Tkachenko an
  fg.satz.kohomologisch; Julia-Gleichung mit Originalquelle an it.lem.julia; Krull/Webster liegen bei L3 (br.frage.bohr_mollerup).

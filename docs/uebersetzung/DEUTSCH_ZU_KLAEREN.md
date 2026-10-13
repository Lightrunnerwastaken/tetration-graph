# Unstimmigkeiten im deutschen Text, bei der Übersetzung aufgefallen

Die Übersetzer sollen deutsche Texte nicht stillschweigend glätten; was ihnen widersprüchlich vorkam, steht hier. Stand 12.10.2026 (Sammlung),
Zweitprüfung 13.10.2026: jeder Punkt fachlich nachgeprüft (eigene Rechnung, wo nötig mit mpmath/numpy). Klare Fehler sind im deutschen Text korrigiert
und in der englischen Fassung nachgezogen (Abschnitt „Erledigt“). Status wurde **nicht** geändert; dazu und zu allem, was ohne Quelle oder Rückfrage
nicht entscheidbar ist, stehen unten Empfehlungen an den Autor.

## Offen — Statusfragen (Empfehlung, nach Wichtigkeit)

- `kn.konstr.pentation` (Status „bewiesen“): λ_b > 1 und Eindeutigkeit des Fixpunkts in (−2, −1) sind nur gemessen. Eigene Überlegung: λ_b ≥ 1 folgt am
  *kleinsten* Fixpunkt aus dem Zwischenwertsatz (T_b − x läuft von −∞ nach 1); λ_b > 1 und Eindeutigkeit folgen, wenn T_b auf (−2, −1) konkav ist — das ist
  gleichwertig zur Log-Konkavität von T_b auf (−1, 0), und die bricht nahe 0 für b > b₇ (`gr.gegen.tpid7_kneser`). **Empfehlung:** Status „bedingt“
  (Bedingung: Log-Konkavität auf (−1, 0) bzw. λ_b > 1), oder den Satz auf die gemessenen Basen beschränken. Das erbt `kn.satz.hexation_struktur`
  (Status „bewiesen“, benutzt die Monotonie von P_b) und mittelbar `kn.konstr.hexation_regulaer`.
- `bw.zert.8r1` (Status „bewiesen“): Prüftext und Bewertung nennen drei offene Punkte (Mikrosegment-Induktion, Hilbert-Röhren-Supremum, Γ-Steigung).
  Die Annahme „korrekt implementiert“ allein wäre bei Rechnerbeweisen üblich, die offenen Punkte nicht. **Empfehlung:** „bedingt“, bis die drei Punkte zu sind.
- `kl.satz.newton_reihe` (Status „bewiesen“): Der Beweis stützte sich auf eine *ganze* Poincaré-Funktion am anziehenden Fixpunkt — die ist dort nicht ganz
  (korrigiert, siehe unten). Gebraucht wird nur Konvergenzradius der Ψ-Reihe > |s₀|; für b = √2 nachgerechnet (Radius ≈ |s₀|λ⁻² = 1.316, |s₀| = 0.632,
  Σ|ψ_k s₀^k| = 1). **Empfehlung:** allgemeines Argument ergänzen (T regulär holomorph auf Re h > −2 und 2πi/ln λ-periodisch ⇒ Ψ holomorph auf |w| < |s₀|λ⁻²),
  sonst Status „bedingt“.
- `kl.satz.carlson` (Status „bewiesen“): Die Eindeutigkeit ist bewiesen; dass die reguläre Tetration für 1 < b ≤ η auf Re z ≥ 0 beschränkt ist, ist gemessen
  (der Kern sagt „gemessen“, f5 fragt danach). **Empfehlung:** Status bleibt; Formulierung „ist nach Messung so eine Lösung“.
- `kn.konstr.carleman` (Status „bewiesen“): Lösbarkeit des Systems für jedes n nur beobachtet. **Empfehlung:** als Definition belassen, aber „Seine Lösung
  (für alle geprüften n eindeutig) definiert slog_n“.
- `gr.satz.superwurzel_dritte` (Status „bewiesen“, Konvergenz der Iteration „nur beobachtet“): Die Konvergenz ist elementar beweisbar — g(y) = exp W(W(y ln Y))
  ist auf [1, ∞) wachsend, g(1) > 1, also steigen die Iterierten ab y = 1 monoton und bleiben unter dem (eindeutigen) Fixpunkt. **Empfehlung:** diese zwei
  Zeilen als Beweis ergänzen und „nur beobachtet“ streichen.
- `gr.satz.wurzel_eta` (Status „bewiesen“, Koeffizient −1.6166 „gemessen“): geschlossen ist er −(5/6)·√(2e/η) = −1.61658… (aus φ(e(1+v)) =
  (1/e)(1 − v²/2 + 5v³/6 + …)), passend zur Messung (−1.61655 bei η − a = 10⁻¹⁰). **Empfehlung:** einsetzen; f5 ist damit beantwortet.
- `fp.satz.messbare_dynamik` (bewiesen_literatur): „nicht ergodisch“ steht nur aus dem Gedächtnis. Meine Erinnerung deckt sich (Rees 1986: nicht rekurrent;
  Lyubich 1987: nicht ergodisch), ist aber ebenfalls keine Quelle. **Empfehlung:** an der zbMATH-Rezension von Lyubich 1987 prüfen.
- `gr.satz.exponentialitaet` (bewiesen; Ganzzahligkeit aus Hardy / van der Hoeven): der Prüftext sagt das schon („(ii) allgemein Literatur“). **Empfehlung:** belassen.
- `kn.satz.hexation_struktur`: gemessene Werte s_1.5 usw. sind als „Gemessen“ markiert — unschädlich; die Statusfrage hängt an `kn.konstr.pentation` (oben).

## Offen — Rückfragen an den Autor

- `it.satz.vertauschbar`: Die Voraussetzung „f^t reell auf der reellen Achse“ ist „benutzt: nein“ (nur Zweigwahl). Ohne sie ist aber auch
  c_t = e^{t(ln λ + 2πik)} eine stetige Gruppe am Fixpunkt holomorpher Iterierter; „die einzige“ gilt also nur mit festem Zweig. Frage: „benutzt: ja“
  setzen (Realität erzwingt c_t > 0) oder „bis auf die Zweigwahl von λ^t“ in die Aussage?
- `fp.satz.p2_raten`, Bewertung: In der Variablen √(e^{−e} − b) ist der Leitkoeffizient jedes Grenzwerts √6·e^{(e−3)/2} ≈ 2.128 (nachgerechnet), nicht
  √(6/e) ≈ 1.486; √(6/e) passt zu einer anderen Normierung (etwa e·z gegen relatives t). Unsere Formel 2√6·e^{−3/2}√t stimmt (1.0931101). Frage: Welche
  Variable benutzt Cho–Kim 2004, Thm. 2.1? Bis dahin „daraus folgt unsere Lücke“ abschwächen.
- `br.satz.boshernitzan`: Kern „(Lösungen von f(x+1) = e^{f(x)})“, Theorie „etwa … am Original zu prüfen“. Frage: Liefert Boshernitzan 1986 wirklich
  eine solche Lösung im Hardy-Feld? Sonst im Kern „(etwa …)“ wie in der Theorie.
- `bw.korr.r20` / `bw.beob.sigma_stern`: dieselbe Zahl 0.058 ist einmal Löserartefakt (Spitzenhöhe), einmal Streifenhöhe aus dem Modenabfall. Die Theorie von
  sigma_stern sagt inzwischen nur noch „stützt …, wenn σ*_∞ > 0 bleibt“. Frage/Empfehlung: σ* mit der exakten Logarithmenleiter (wie in r20) neu messen,
  bevor es `bw.verm.streifen` stützt.
- `gr.satz.frei_und_fest`: Theorie an den Kern angeglichen („endlich viele Ableitungen“). Für C^∞-Tetrationen ginge nach Borel sogar jede Taylorreihe an einer
  nicht ganzen Stelle (die ganze Stelle hat schon „alle Ableitungen“). Frage: Ist die vorsichtige Fassung gewollt?
- Konditionszahl `kn`: „≈ 1.9 je Grad“ (Theorie `schranke`, f4 von `kn.satz.polynom_halbexp_endlich`) und „≈ 2.6 je Grad“ (`β_N` in `kantorovich`)
  messen erkennbar verschiedenes (Konditionszahl an p_N gegen ‖J⁻¹‖ an Knesers q_N). Bitte beide Grössen ausdrücklich benennen.
- Quelle `justus_bc_geteilt`: deutscher Titel „Base Change an einem geteilten Fixpunkt“ — gemischtsprachig; Originaltitel der Notiz?
- `quellen[].stelle`: mehrere Einträge sind Kommentarsätze statt Fundstellen; „Levenstein 31.03.2018 Grenzfunktion für Kneser 2 gegen e“ ist mehrdeutig.
  Kein Fehler, aber für die Ansicht unpraktisch; bei Gelegenheit in `notiz` verschieben.
- Bitte gegenlesen (Deutung): Strang `kl.gleiches_L` „bleibt er stehen“ habe ich als „läuft der naive Höhenvergleich gegen eine Konstante“ gelesen
  (so `exp.koenigs_basiswechsel`, Prüfung 5); Szekeres: auch f5 von `it.satz.szekeres` („wandte ein, Knesers Lösung …“) auf 1962 gesetzt.

## Erledigt (2026-10-13)

**fg.** `fg.satz.allgemeine_loesung`: Titel „Alle stetigen, streng wachsenden Lösungen“ (wie Bereich). `fg.satz.tpid4`: kein Widerspruch — der Satz
beweist die Vermutung (Fourier-Argument nachgerechnet).

**it.** Szekeres: einheitlich 1962 (`it.satz.szekeres_monoton` Beweis und Theorie, f5 von `it.satz.szekeres`, Lehrbuch Kap. 0); die Quelle selbst
datiert J. Austral. Math. Soc. 2 auf 1962 (Band 1961/62), Schlüssel `szekeres1961` bleibt. `it.satz.cremer`: kein Fehler — JFM 53 (Jahrgang 1927) und
Math. Ann. 98 (1928) stehen beide in der Quelle. `it.satz.parabolische_rate`: Summationsindex j statt k (Kollision mit dem Exponenten k; Rechnung
nachgeprüft, auch das η-Beispiel). `it.satz.vertauschbar`: „stetig in t“ in den Bereich (Zweigfrage oben offen).

**fp.** `fp.beob.i_turm`: f3 „b = i“; Beweis nennt jetzt fp.satz.shell_thron (t = λ, te^{−t} = ln i, |t| = 0.8915; numerisch bestätigt) statt „1 im Becken
nur numerisch“ — Status bleibt. `fp.beob.regionen`: Perioden 5–12 zusammen 42 198 ergänzt (Summe 270 400, aus dem Experimentergebnis).
`fp.beob.partner_keil`: Keilspitze auf der negativen reellen Achse nahe 0 statt „bei e^{−e}“ (nachgerechnet: alle Basen haben Re b ≤ −0.0064), Index
einheitlich m. `fp.verm.knospen`: Kennung historisch, bleibt. `fp.satz.rechteckwelle`, `fp.satz.kreismittel`: kein Widerspruch — die Voraussetzungen stehen
als „benutzt: nein“. `fp.satz.siegel_mittel`: kein Widerspruch (Satz für beschränkte Scheiben, Skizze für die wirklichen; Status „skizze“ passt).

**br.** `br.satz.ordinal_tetration` einfach: „ω hoch ε₀“. `br.frage.transexp`: kein „\num“ im Quelltext (Zeilenumbruch + „um“ im JSON-Export).
`br.satz.nicht_definierbar`, `br.satz.young_kreis`: „benutzt: nein“, konsistent. `br.satz.halbexp_komplexitaet`: Bereich (Definition f∘f ≈ 2^n),
Aussage (technische Bedingung 2^{o(n)}) und Theorie sagen Verträgliches. `br.satz.shishikura`: `\text{abst}` → `\text{rep}` ist erlaubte Übersetzung von `\text`-Prosa.

**ar.** `ar.lem.nullstellen` f2: Holomorphie auf D̄(1, 2.9) (jeder Radius > 1.2 ginge, andere Konstanten). `ar.kor.benachbart`: x = 0 ist keine Ausnahme
(T(1) = e transzendent; der Beweis deckt x = 0), Aussage und Theorie „x ≠ −1“. `ar.satz.selbstpotenz_kongruenz` einfach: nur teilerfremde Reste
(mod 4 kommt 2 nie als x^x vor). `ar.verm.gelfond_turm`: Bereich „ω ≠ 0 algebraisch“, Theorie mit ¹z. `ar.satz.unendliche_tuerme_gs`: Bereich ab e^{−e}
(Gelfond–Schneider-Argument trägt dort), Theorie ebenso. `ar.frage.ef_fortsetzung`: einfach ist zulässige Vereinfachung. `ar.satz.ssc_gelfond`: die
Implikation ist bewiesen (anders als `ar.satz.schanuel_orbit`, das die Folgerung selbst behauptet) — Status passt.

**gr.** einfach: `gr.satz.grenzfaelle` (Höhe −2), `gr.satz.gestalt` („wachsende“), `gr.satz.integral_produkt` („zwischen ganzen Zahlen“).
Theorie an Kern: `gr.satz.stellen` (ganze a ≥ 2), `gr.satz.logstern` (Bereich y > 0 wie die Theorie — die Identität gilt dort; worum „höchstens 1“),
`gr.satz.superwurzel` (Fall x = 1), `gr.satz.frei_und_fest` („endlich viele“), `gr.satz.turmintegral_periodisch` („genau einmal“ als Voraussetzung;
„linear“ = x^s ist die übliche lineare Tetration 1 + s auf [−1, 0], nicht vertauscht). `gr.beob.wendepunkte`: −1.94 (Skript e12). `gr.beob.kneser_pde_wackeln`:
Titel „±3·10⁻⁴“ (halbe Spanne 5.7·10⁻⁴). `gr.satz.basis_ganz`: „der Kneser-Familie“. `gr.verm.tpid7`: Theorie — der Wendepunkt wandert mit der Basis nach
**links** (b = 10: −0.997), nicht „zu nah an 0“; „folgt aus der Vermutung kl.verm.vollmonoton“; „höchstens bis b₇“. `gr.gegen.tpid7_kneser` Theorie:
„bis b₇ nur gemessen“. `gr.satz.tauschregel` f3: Exponentialfunktion. Kein Fehler: `gr.satz.young` („benutzt: nein“), `gr.satz.stufenleiter`
(abgeschnitten gegen gerundet), `gr.satz.kehrwert_turmsumme` (je nach Basis das vierte oder fünfte Glied).

**kl.** `kl.satz.newton_reihe`: Beweis und f2 ohne „ganze Poincaré-Funktion“ (Statusfrage oben). `kl.satz.invariante_daempfung`: Grenzfall |λ| → 1 als
Deutung (wie Lücke und Theorie). `kl.beob.zwei_abel`: 170 Stellen, 3.6·10⁻²⁵ (Experiment e03), Titel „periodische Funktion“.
`kl.beob.stetig_eta`: fünf Stellen. `kl.satz.koenigs_basiswechsel`: Konjugation 10⁻³⁰ (Experiment 5.9·10⁻³¹), Höhe 10⁻²⁷ — verschiedene Grössen.
`kl.beob.eta_konstanten`: Tiefe 3000 (Fatou-Koordinate) und 20 000 Turmstufen (C_T, Summen) — beide richtig, im Prüftext unterschieden.
`kl.satz.bruchsumme`: Dezimalpunkt. `kl.satz.basiswechsel_gleiches_L`: „benutzt: nein“, konsistent.

**bw.** `bw.prop.saddle`: Konturverschiebung in den Bereich (dann ist „bewiesen“ unstrittig). `bw.satz.kern_existenz`: (K3) = Holomorphie in der Basis
(so `kn.satz.paulsen_komplex`), Verweis ergänzt; Thm. 13 (Existenz) und Thm. 2 (Eindeutigkeit) widersprechen sich nicht. `bw.beob.m2`: Levenstein-Werte als
slog_3 − slog_e ≈ −0.1926 = −m(3) usw.; „vier Stellen“ (signifikant) stimmt. `bw.gescheitert.vorhersage_r8`: Titel „Exponent 2πν“.
`bw.satz.imaginaere_richtung`: Dämpfung je Umlauf 3.5·10⁻¹⁵ (nachgerechnet); Bereich Kern/Theorie umschreiben dieselbe Bedingung (f2).
`bw.satz.paulsen_cowgill_hardy`: Status nach Repo-Konvention (sekundär geprüft); die Glosse ist die Definition. `bw.korr.konstante`: „2.3 %“ von der Korrektur
unberührt, selbst aber zu optimistisch (4–7 %, `bw.lem.randgesetz`). `bw.beob.kernkleid`: bei d ∈ {2, 10} RH-Löser statt komplexer Schritt.
`bw.verm.nirgends_analytisch`: „gleichmässige Konvergenz aller Ableitungen“. `bw.satz.fourier_naht`: „Sprung“. `bw.prop.antisym_c1`: „benutzt: nein“, konsistent.

**kn.** `kn.satz.polynom_eckprofil`: „(i) bewiesen bis auf Lemma K12“ (Status „skizze“ passt). `kn.verm.tpid11` einfach: drei nachgerechnete Stufen
(Abstände 1.6 %, 2.3 %, 2.0 %). `kn.verm.tpid3`: die drei anderen Fälle (Rand, ausserhalb, Polynom p = 2). `kn.satz.walker_konjugation` f4: Konjugation
von b^x zu einer Abbildung mit parabolischem Fixpunkt; zwei Walker-Arbeiten 1991 (Math. Comp., J. Math. Anal. Appl.) stehen beide in den Quellen.
`kn.beob.hoehere_etas`: Theorie nennt „genau bis ηₙ“ als Vermutung (`kn.verm.beschraenkt_unter_2`); einfach „einziger Wert für das fünfte Glied“.
`kn.satz.paulsen_komplex`: Theorie an den Kern (Holomorphie in b nicht belegt). `kn.konstr.hexation_regulaer` f1: „zwei Stufen höher“.
`kn.frage.hexation_kneser` Theorie: das Fork-Verfahren rechnet nie mit θ auf der Achse. `kn.verm.polynom_halbexp` Theorie: bedingt bewiesen, offen ist
`kn.verm.polynom_halbexp_kern`. `kn.satz.halbexp_beschraenkt`: Voraussetzung in die Aussage. `kn.satz.polynom_halbexp_hauptsatz` einfach: „mit wachsendem Grad“.
Kein Fehler: `kn.satz.polynom_approximation_ecke` („sheet rotation“ schon richtig), `kn.satz.kneser_e_zertifiziert` (Streifen im Bereich definiert),
`kn.satz.abel_bruchsumme`, `kn.satz.singulaer_minus2` (je „benutzt: nein“).

**Lehrbuch und README.** Kapitel 0: „bei höchstens 1“. Formelsammlung: „bei a = e kein Sprung“. README (und README.en.md): „9 Gebiete“ stimmt
(11 Dateien, 9 Gebiete; die zehn Zahlen sind Kapitel — „damals je Kapitel“); `formelnEinpassen` nur noch einmal beschrieben; e38 (`e38_hexation.py`)
im Hexations-Absatz, der jetzt eine eigene Überschrift (09.10.2026) hat; dabei gefunden: b₇ stand als 10.2210634581 statt 10.2210634158.
Strukturtexte: `fp.turm` (Verb, reell von e^{−e}), `fp.parameterraum` („kommt der Turm zur Ruhe“), `kl.gleiches_L` (siehe oben),
`exp.polynom_team` („zu Π_N f_K“), `kn.hoehere_stufen` („und, anders als Knesers Tetration, … verzweigt“).

**Quellen, Experimente.** `forum_tid1177`: „Nuninho“. `euler1778`: Dezimalpunkt. `forum_tid213`: sagt schon „1941 … besprochen“.
„Abschnitt Super-root“ ist der Abschnittsname der englischen Wikipedia. `exp.ableitungen_vorzeichen`: „20 Stellen“ = 20 Gitterpunkte
(Prüfungsname); die englische Fassung sagt schon „points“, die deutsche Zusammenfassung erzeugt das Skript.

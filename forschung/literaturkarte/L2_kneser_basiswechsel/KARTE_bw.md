# Literaturkarte bw — Basiswechsel (Tetrationen verschiedener Basen ineinander umrechnen)

Strang L2, Stand 12.10.2026. Die eigene Serie (Papers I–VI, Synthese, Zertifikat 8r1, imaginäre Richtung) ist in 06_basiswechsel.yaml
ausführlich; hier geht es um die Arbeiten **anderer**. Suchwege: Tetration-Forum (Suche „base change“, „change of base“, Thema 1189 und seine
Querverweise, Thema 262 Seiten 2–5), zbMATH Open, Crossref, Korean J. Math. (Volltext Paulsen 2016), MathOverflow/Math.SE (über Forumsverweise).
Prüfstand P wie in KARTE_kn.md (O Original, S sekundär, N nicht geprüft).

## 1. Die Basiswechsel-Formel und ihre Konstante

| Wer, wann | Wo | Ergebnis | Methode | Beweisstand | P | Bei uns |
|---|---|---|---|---|---|---|
| J. D. Fox 12.08.2007 | Forum tid 14 „Change of base formula for Tetration“ | ^x a = lim log_aⁿ(^(n+x+μ_b(a)) b), additive Konstante μ mit μ_a(c) = μ_a(b) + μ_b(c) | Grenzwert iterierter Logarithmen | Heuristik; μ als exakt angenommen | O | forum_fox2007; bw.def.F_bd, bw.prop.hub |
| G. Helms 16.08.2007 | Forum tid 26 „Gottfried Helms' change of base formula“ | Basiswechsel zwischen **regulären** Tetrationen unter η (Basen t^{1/t}, 1 < t < e) über Eigenzerlegung der Carleman-Matrix; Beispiel √2 → 2.5^{1/2.5}, y₂ ≈ 2.342 | Eigenzerlegung abgeschnittener Matrizen | numerisch, ungenau | O | — → L1 (kl.satz.koenigs_basiswechsel) |
| S. Levenstein 18.–19.02.2009 | Forum tid 236 | slog₂ − slog_e → ≈ 1.1282, slog₃ − slog_e ≈ −0.1926, slog₁₀ − slog_e ≈ −1.1364; „wackelt ein wenig“ | Tabellenrechnung | numerisch | O | forum_levenstein2009; bw.beob.m2 |
| H. Trappmann, S. Levenstein Aug. 2009 | Forum tid 333 „Real and complex behaviour of the base change function“ | log_aⁿ∘exp_bⁿ hat Verzweigungspunkte bei log_bⁿ{a, …}; nur reell C^∞ | Analyse der Verzweigungen | Skizze | O | forum_trappmann2009 |
| BenStandeven 30.04.–01.05.2009 | Forum tid 262, Beiträge 22–27 | Fixpunktstörung für Basis η+ε; Engpass: slog_η(e−√ε) ≈ 2e/√ε; Turm (η+ε)↑↑(C/√ε) → e | Störungsrechnung | Heuristik | O | — (neu bei bw.lem.randgesetz) |
| S. Levenstein 20.01.2011 (dazu tid 825, 2014; Math.SE 2010406) | Forum tid 578 „Nowhere analytic superexponential convergence“ | Vermutung: Basiswechsel-Tetrationen (und Walkers Lösung) sind C^∞, aber nirgends analytisch; Konvergenzradius springt von ≈ 0.46 auf ≈ 0.035 | Taylorkoeffizienten bis Millionen | Vermutung, ungelöst | O | forum_levenstein2011; bw.verm.nirgends_analytisch |
| Gottfried 03.10.2012 | Forum tid 750; Math.SE 198757 „A tetration limit for base a = e^(1/e)“ | Grenzwertfrage zum Basiswechsel an η | — | Frage | S | — → L1 (Basis η) |
| tommy1729 2012 | Forum tid 755 „base change inconsistant ?“ | Zweifel an der Ableitung von b↑↑x über den Basiswechsel | — | Diskussion | N | — |
| J. D. Nixon, S. Levenstein 2014 | Forum tid 937 „Fractional calculus and tetration“ | Basiswechsel über η hinaus: alle Ableitungen auf der reellen Achse konvergieren, wachsen aber zu schnell für Analytizität; im Komplexen nicht definiert | Ableitungsabschätzung | Heuristik | S | — (bei bw.verm.nirgends_analytisch ergänzen) |
| S. Levenstein 11. und 15.06.2015 | Forum tid 262, Beiträge 40–41 | **Randgesetz**: sexp_{η+1/n}[π√(2ηn/e) + slog_e(n/e − 1) + C] − n → 0 (Engpass-Durchgang wie an der Spitze der Mandelbrot-Menge); C ≈ −2 − slog_e(388.7874/e − 1) | Engpass-Rechnung, Analogie zu π in der Mandelbrot-Menge | Vermutung mit Heuristik | O | — (neu bei bw.lem.randgesetz) |
| W. Paulsen 2016 | Finding the natural solution to f(f(x)) = exp(x), Korean J. Math. 24, 81–106 | Stetigkeit seiner Lösung bei a → η (Korollar 1); im Beweis von Lemma 2 die Engpass-Zählung π√2/√(η_a f''(x_a)) für Iterationen nahe einer Fast-Berührung | Störung + Engpass-Riemann-Summe | bewiesen | O (S. 102–105) | paulsen2016 (bisher nur in kn) |
| D. Taranovsky (MO 2017), S. Levenstein 28.–31.03.2018 | Forum tid 1189 | Ordnen sich gebrochene exp-Iterierte zu verschiedenen Basen? Für Kneser nein (Welle), für Walker vermutlich ja; Levenstein verweist auf die kommende Arbeit von Paulsen–Cowgill | Wellenbild slog_e(tet₂(x)) − x | Vermutung | O | forum_tid1189; bw.frage.halbiter_vergleich |
| W. Paulsen, S. Cowgill 2021 | Incompatible growth rates of the tetrations, Indian J. Math. 63(1), 127–141 (zbMATH 1470.39046) | Kneser-Tetrationen zu zwei verschiedenen Basen liegen nicht in einem gemeinsamen Hardy-Körper; auch ihre gebrochenen Iterierten sind unverträglich | Hardy-Körper, Eindeutigkeit der Tetration | bewiesen | S (zbMATH-Zusammenfassung; Text nicht zugänglich) | neu: bw.satz.paulsen_cowgill_hardy |
| J. D. Nixon 17.03.2021 | Forum tid 1308 „A different approach to the base-change method“ | Basiswechsel von η nach e über Mellin-Darstellung der gebrochenen Iterierten zur Basis η; reelle Wahl von θ offen | Mellin-Transformation, Schröder an L | Skizze, offen | O | — |
| P. L. Walker 1991 | Math. Comp. 57, 723–733 | C^∞-Lösung als Grenzwert dexp^{−n}∘expⁿ — ein Basiswechsel von der parabolischen eˣ−1 zu exp | Grenzwert | C^∞ bewiesen | S | kn.konstr.walker |
| W. Paulsen 2026 | Analyzing the function z↑↑a for fractional a, Adv. Comput. Math. 52(2), Nr. 25 | Tetration fester Höhe als Funktion der Basis | — | — | S (nur Metadaten) | — (Lücke A) |

## 2. Wer hat welche Frage verschieden angegangen

- **Ist der Basiswechsel eine Konstante?** Fox 2007 nahm μ als exakt an; Levenstein 2009 sah das Wackeln in Tabellen; Levenstein 2018 zeichnete die
  1-periodische Grenzfunktion (Kneser 2 gegen e); Paulsen–Cowgill 2021 beweisen die Hardy-Körper-Unverträglichkeit; unsere Papers I–IV
  fassen es als Phasengesetz mit Kreishomöomorphismus (bw.satz.phasengesetz, bw.satz.homoeo).
- **Ist die Basiswechsel-Tetration analytisch?** Trappmann–Levenstein 2009 (Verzweigungspunkte), Levenstein 2011/2014 (Konvergenzradius
  bricht ein: nirgends analytisch, vermutet), Nixon–Levenstein 2014 (Ableitungen wachsen zu schnell); unsere Streifenvermutung bw.verm.streifen
  sagt das Gegenteil für die Phase — der Widerspruch steht offen in bw.frage.streifen_gevrey.
- **Verhalten nahe η:** BenStandeven 2009 (Engpass 2e/√ε), Levenstein 2015 (Konstante π√(2η/e) mit log-artiger Korrektur slog_e(n/e−1)),
  Paulsen 2016 (Engpass-Zählung im Beweis), Paper IV Lemma 8.2 (bw.lem.randgesetz) — alle dasselbe Sattel-Knoten-Gesetz.
- **Basiswechsel im Komplexen:** Trappmann–Levenstein 2009 (Verzweigungen), Nixon 2021 (über Mellin), unsere imaginäre Richtung
  (bw.satz.imaginaere_richtung: Grenzwert A_b(L_d)) — der Grenzwert in imaginärer Richtung wurde im Forum nicht gefunden.

## 3. Lückenliste bw (Priorität A = zuerst)

- **A** Paulsen–Cowgill 2021 als Knoten; Bezug zu bw.satz.phasengesetz und bw.frage.halbiter_vergleich (beantwortet sie für Kneser?) — Text
  besorgen (Indian J. Math., nicht frei); Leitung fragen.
- **A** bw.lem.randgesetz: Neuheit von „folklore“ auf „bekannt“ (Levenstein 2015, Paulsen 2016, BenStandeven 2009). → in diesem Durchgang.
- **A** Paulsen 2026 lesen (Basisabhängigkeit, bw.ann.basis_diffbar).
- **B (erledigt 12.10.)** Helms 2007 an L1 gemeldet (BRETT); Quelle `forum_tid26` liegt in lit_L2.yaml.
- **B (erledigt 12.10.)** Nixon 2021 (tid 1308) → bw.frage.nixon_eta_basiswechsel; Nixon–Levenstein 2014 (tid 937) als Quelle an bw.verm.nirgends_analytisch.
- **C (nicht erreichbar)** Math.SE 2010406 und 198757: Abruf von math.stackexchange.com gesperrt.

## 4. Stand und Restlücken (12.10.2026, nach Block 3)

**Stand:** Alle gefundenen Arbeiten anderer zum Basiswechsel sind kartiert und mit unseren Knoten verbunden: Fox 2007 (Formel), Levenstein 2009/2018 (Welle),
Trappmann–Levenstein 2009 und Levenstein 2011/2014 (Analytizität), BenStandeven 2009, Levenstein 2015 und Paulsen 2016 (Randgesetz bei η → bekannt),
Paulsen–Cowgill 2021 (Hardy-Körper, bw.satz.paulsen_cowgill_hardy), Nixon 2021 (Basiswechsel von η aus, bw.frage.nixon_eta_basiswechsel). Neu bewertet:
bw.lem.randgesetz (bekannt), bw.satz.imaginaere_richtung (vermutlich neu, Recherche dokumentiert).

**Restlücken:** Volltexte Paulsen–Cowgill 2021 und Paulsen 2026 (Bezahlschranke, bleiben sekundaer) — entscheiden, ob bw.frage.halbiter_vergleich für Kneser beantwortet ist
und ob die Basisabhängigkeit (bw.ann.basis_diffbar) belegt werden kann; Math.SE 2010406/198757 nicht erreichbar; die übrigen bw-Aussagen sind eigene Forschung
(Papers I–VI) ohne Gegenstück in der Literatur.

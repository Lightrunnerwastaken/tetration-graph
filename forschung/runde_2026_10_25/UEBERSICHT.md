# Forschungsrunde 25.10.2026 — Übersicht

Geplant waren 70 Gebiete mit höchstens 16 gleichzeitig laufenden Agenten. Der Lauf startete um 13:20 und wurde auf Wunsch um 14:00 angehalten.
**18 Gebiete sind mit Bericht fertig** (01–16, 23, 29). 14 Gebiete wurden mitten in der Arbeit abgebrochen und haben keinen Bericht (17–22, 24–28, 30–32;
ihre Teilskripte liegen in den Ordnern). Die Gebiete 33–70 sind nicht gestartet. Die Liste aller 70 Aufträge steht im Workflow-Skript der Sitzung.
Die Befunde sind von den Agenten eingeordnet; ich habe sie nicht nachgerechnet. Vor dem Eintragen in den Graphen ist eine Gegenprobe nötig.

## Stärkste Befunde

| # | Befund | Status laut Bericht | Ordner |
|---|---|---|---|
| 1 | **Zweigwahl der θ-Singularität mit Intervallarithmetik bewiesen** (Basis e und 2): es genügt P ∈ σ(ℍ), von den Motordaten wird nur θ(z₀) auf 0.027 verlangt (gemessen 1e-59). Damit werden kn.satz.theta_singularitaet und kn.satz.theta_koeffizienten von „bedingt“ zu „bewiesen“ (unter dieser schwachen Datenannahme). | bewiesen | 03_zweigwahl_zertifikat |
| 2 | **Natürliche Grenze bewiesen:** Für jede Superfunktion an einem anziehenden Fixpunkt ist der Rand des Becken-Urbilds eine natürliche Grenze der θ-Welle (nur Koenigs, Identitätssatz, Rechtsverschiebung). Löst die Vermutung zur Kneser-artigen Pentation aus Runde 18. | bewiesen | 08_pentation_natuerliche_grenze |
| 3 | **Monodromie des regulären Erzeugers vollständig:** Blattadressen in ℤ^(∞); die Monodromiegruppe ist die iterierte Monodromiegruppe IMG(exp), ein iteriertes Kranzprodukt ℤ≀ℤ≀…, nicht frei. | bewiesen | 15_erzeuger_monodromie |
| 4 | **Explizite Phase der θ-Koeffizienten:** arg aₙ = π + arg(1+θ′(σₙ)) + arctan(π/(2(ln n − ln(A/2π)))) + O(1/(ln n ln ln n)), kompakt aₙ ≈ −(1+θ′(σₙ))/(2πn²) mit σₙ ≈ 1/(2πn) (weil 1+θ′ dort etwa 2πn/ln n ist, passt das zu |aₙ| ~ 1/(n ln n)). Beantwortet den letzten offenen Punkt von kn.frage.theta_exponent. | hergeleitet mit Lücke + numerisch | 02_theta_argument |
| 5 | **Rang-5-Basiswechsel-Phase ist singulär (Cantor-Typ)**, nicht Weierstrass-artig: streng wachsend, Ableitung 0 fast überall, Hölder-Exponent trotzdem überall 1 (Kakutani/Riesz-Produkt). **Widerspricht** der Vermutung „nirgends differenzierbar, Weierstrass-artig“ aus Runde 18 (bw.verm…). | hergeleitet mit Lücke | 10_basiswechsel_hoelder |
| 6 | **Gesetz für gerade Ränge hergeleitet** aus einem lokal bewiesenen Kern-Lemma H_s⁻¹(V) = z_s + Φ(λ_{s−2}^V); bei Rang 10 neu bestätigt: α₁₀ gemessen 0.8360, vorhergesagt 0.83607. Offen nur die Zweigwahl (vgl. Befund 1). | hergeleitet mit Lücke + numerisch | 09_gerade_raenge_gesetz |
| 7 | **Riemannfläche von sexp:** Verzweigungsorte aller Nebenblätter = {z : sexp(z+M) = sexp(M−n), n ≥ 2}; dicht in ℂ, abzählbar, für jede feste Tiefe nirgends dicht. | bewiesen (Dichtheit mit Lücke) | 06_riemannflaeche_sexp |
| 8 | **Riemannfläche von slog:** links von Re w = Re L ein 2πi-periodisches Blatt mit Log-Gitter L+2πiℤ, L̄+2πiℤ; rechts Knesers Blatt, auf dem diese Punkte regulär sind (z. B. sexp(3.0339794527621306885 + 0.6260388922275826520i) = L); exakte Monodromie an der Naht. | bewiesen + numerisch | 07_riemannflaeche_slog |
| 9 | **Q* scharf eingegrenzt:** Q* gilt für jede auf ℍ holomorphe Tetration, deren Kneser-Abweichung p_T eindeutig in die obere Halbebene reicht und o(e^{2π Im z}) wächst; Gegenbeispiele brauchen Anteile e^{−2πinz} oder Mehrdeutigkeit. Schwarz-Starrheit: halbebenen-erhaltende Verbiegungen sind trivial. | 8 Aussagen bewiesen | 01_halbebene_eindeutigkeit |
| 10 | **Momentenfamilie:** gleichwertig zu „Cauchy-Transformierte der Werteverteilung = 1/(ζ−c) ausserhalb der Kurve“; damit ist jede Tetration mit doppelpunktfreier θ_T-Kurve, die die Familie erfüllt, Kneser; ausdrückliche Umgebung über ein Sehnenkriterium. | bewiesen/hergeleitet | 04_momentenfamilie |
| 11 | **Ausgewogene Tetration = linke Pentation = Iteration von m·e^m** (Umkehrung von Lambert W); alle linken Ränge ab 5 parabolisch bei 1 mit geschlossenem iterativem Residuum ((9−2k)/(2x)), daher eindeutige natürliche Fortsetzung. | bewiesen | 14_klammerungen |
| 12 | **Fluchtmuster von eᶻ in Kneser-Koordinaten doppelt periodisch:** Gitter ℤ + (2πi/L)ℤ, Eintrittsindex N(z+ω) = N(z)+4 mit ω = 2πi/L − 4; für die reguläre Superfunktion exakt bewiesen. | bewiesen (regulär) / hergeleitet (Kneser) | 16_julia_menge_sexp |
| 13 | **Reguläre und Kneser-Hexation stimmen an η₅ nicht überein:** sie unterscheiden sich um eine 1-periodische Höhenverschiebung mit Amplitude ~1e-5 (fast nur erste Harmonische). | numerisch | 12_hoehere_raenge_basis |
| 14 | **Iteratives Residuum ι_r der Berührstellen** als neue Kennzahl der Rangleiter (r = 3…8 auf ≥ 17 Stellen) mit bewiesener verfeinerter Rate; Vermutung 2 − η_r ~ 2.20/(r ln r). | numerisch + bewiesen | 11_eta_n_asymptotik |
| 15 | **Transzendenz:** Σ 1/(Turmwerte) ist transzendent (Kettenbruch explizit; im Graphen bisher „offen“); T′(0), T′(1) nicht beide algebraisch; bei ganzer Basis höchstens ein T′(n) algebraisch; 47 Kneser-Werte mit 21 485 Kettenbruchgliedern ohne Struktur. | bewiesen + numerisch | 29_transzendenz |
| 16 | **Basisableitung:** festes Vorzeichenmuster für b von 1.000001 bis 10⁴ (steigt auf (−2,−1), sinkt auf (−1,0), steigt ab 0); Kurven verschiedener Basen schneiden sich nur bei −1 und 0; alles hängt an einer „Elastizität“ auf (−1,0); Kanalgesetz (b−η)^{−3/2} auf 6e-5. | bewiesen + numerisch | 23_basisableitung |
| 17 | **Komplexe Basen:** nahe dem Shell–Thron-Rand wandert eine zusätzliche Singularität der Welle in die obere Halbebene, das Koeffizientengesetz bricht dort. | numerisch | 05_komplexe_basen_welle |
| 18 | **Rangkoordinate:** Φ = 3.5 legt die Halbrang-Leiter weit oben fest, aber nicht H₃.₅ selbst (zwei verschiedene Familien bei √2). | numerisch, negativ | 13_rangkoordinate |

## Widersprüche und Zusammenhänge
- **10 gegen Runde 18 (01):** „Weierstrass-artig, nirgends differenzierbar“ (Vermutung 18) gegen „singulär, Ableitung 0 fast überall, Hölder 1“ (25/10). Beide hergeleitet mit Lücke — vor dem Eintragen klären; der Graph-Eintrag zur Vermutung muss angepasst werden.
- **03 stützt 09 und kn.satz.theta_***: die zertifizierte Zweigwahl schliesst die offene Bedingung in mehreren Einträgen.
- **06, 07, 15** beschreiben dieselbe Struktur (Riemannflächen von sexp, slog und dem Erzeuger) von drei Seiten — beim Eintragen zusammen lesen.
- **12 und 13** setzen die Beobachtung „kleine 1-periodische Wellen zwischen Konstruktionen“ auf höheren Rängen fort.

## Neue Gebiete
- Monodromie/iterierte Monodromiegruppe von exp im Zusammenhang mit Tetration (15).
- Dynamik von eᶻ in Kneser-Koordinaten: doppelt periodisches Fluchtmuster (16).
- Singuläre (Cantor-artige) Phasenfunktionen beim Basiswechsel höherer Ränge (10).

## Offene Fragen für eine nächste Runde
1. Widerspruch 10 ↔ Runde-18-Vermutung zur Rang-5-Welle auflösen.
2. Dichtheitsbeweis der Nebenblatt-Verzweigungen von sexp vervollständigen (06).
3. Q* für mehrdeutige oder schnell wachsende p_T entscheiden (01).
4. Phasenformel für arg aₙ streng machen (02).
5. Fluchtmuster-Formel für Kneser (nicht nur regulär) beweisen (16).
6. Vermutung 2 − η_r ~ 2.20/(r ln r) an η₉, η₁₀ prüfen (11).
7. Ursache der Hexations-Welle an η₅ erklären (12).
8. Die 14 abgebrochenen und 38 nicht gestarteten Gebiete (17–70) nachholen.

# Eintrag Bericht 02 in den Graphen

Dateien: `daten/aussagen/r18_02_kn.yaml`, `daten/theorie/r18_02_kn.yaml`.

## Angelegte ids
- kn.satz.theta6_singularitaet | θ der Hexation singulär bei x* = ihex(p) + ℤ, Wurzelknick α = 0.4218 (B1–B3)
- kn.satz.theta6_koeffizienten | aₙ ≈ C₆ e^(−2πinx*) n^(−1−α), Modell mit Rest O(1/n²) (B4, B5)
- kn.verm.theta_gerade_raenge | Gesetz für gerade Ränge, Ort H_r⁻¹(p_{r−2}), Exponent ln λ_{r−4}/ln λ_{r−2}; Rang 8 (B6, B7)
- kn.beob.theta5_wachsende_koeffizienten | Rang 5 am anziehenden L₄: Koeffizienten wachsen, θ reicht nicht bis ℝ (B8, B9)
- kn.verm.theta5_natuerliche_grenze | Beckenrand von L₄ als natürliche Grenze (B9)
- kn.beob.pentationswelle_singularitaet | Welle regulär ↔ Kneser-artige Pentation, Abfall e^(−2π·1.6902n)/n, Rangskalierung (B10, B11)
- kn.beob.theta_basis_quasiuniversal | |aₙ(b)| für kleine n fast basisunabhängig, asymptotisch nicht (B12, B13)

Nicht eigens angelegt:
- `kn.beob.hexation_zwei_fixpunkte` (B1): Beide Fixpunkte der Pentation mit denselben Werten stehen schon in `kn.satz.kein_kneser_anziehend`; B1 ist im Bereich von kn.satz.theta6_singularitaet vermerkt.
- Exponent 0.4218 und Ecke der Kurve sind qualitativ schon in `kn.beob.hexation_verzweigt` (x_c = x* + 1); neu sind Ort als ihex(p), Konstanten, Koeffizientengesetz. Kante `stuetzt` gesetzt.
- Ergänzung zu `kn.satz.theta_singularitaet` (f5) nicht in die fremde Datei geschrieben: Vorschlag, dort in f4 auf kn.satz.theta6_singularitaet zu verweisen.

## Vorschlag Strang
Bestehenden Strang zu Knesers θ-Welle (θ_singularitaet → θ_koeffizienten) fortsetzen, oder neuer Strang **„Knesers Welle auf allen Rängen“** — verlauf: Rang 4 (kn.satz.theta_singularitaet, kn.satz.theta_koeffizienten) → Rang 6 (kn.satz.theta6_singularitaet, kn.satz.theta6_koeffizienten) → Gesetz gerade Ränge (kn.verm.theta_gerade_raenge) → Rang 5 ohne reelle Welle (kn.beob.theta5_wachsende_koeffizienten) → Welle zwischen den Pentationen (kn.beob.pentationswelle_singularitaet).

## Vorschlag Ideenkarte
Unteridee bei der Kneser-θ-Idee bzw. bei den höheren Rängen: **„Ausnahmewert bestimmt die Singularität“** — kurz: Die Welle ist dort singulär, wo die Kneser-Lösung einen Wert trifft, den die reguläre Lösung auslässt; Art des Knicks aus dem Verhältnis zweier Multiplikator-Logarithmen.

# Eintrag Bericht 12 (Höhere Ränge als Funktion der Basis) — 25.10.2026

Dateien: `daten/aussagen/r25_12_kn.yaml`, `daten/theorie/r25_12_kn.yaml`. `python baue.py --nur-pruefen`: 0 Fehler, keine Warnung zu diesen ids.

## Angelegte ids

| id | Titel kurz | Status |
|---|---|---|
| kn.satz.ganze_negative_hoehen | H_r(b,−k) = 1 − k für 1 ≤ k ≤ r − 3 (B1) | bewiesen (Voraussetzung Injektivität), Neuheit bekannt |
| kn.satz.negative_struktur | gerade Ränge: Pol σ_r und Fixpunkt p_r in (−(r−2), −(r−3)); ungerade: kein negativer Fixpunkt (B2) | bedingt (Monotonie) |
| kn.verm.singulaere_basen_ungerade | nur η₃, η₅, η₇, … sind singuläre Basen; η₄, η₆, … nur Schwellen (B3) | vermutung |
| kn.beob.pentation_analytisch_eta4 | H₅ an η₄ analytisch in b, Cauchy-Koeffizienten auf 10⁻¹⁹ (B4) | numerisch |
| kn.beob.hexation_sprung_eta5 | Hexation springt an η₅: R(x) = L(x + θ(x)), θ 1-periodisch ~10⁻⁵ (B5) | numerisch, Deutung offen |
| kn.beob.landkarte_p4_p6 | p₄, σ₆, p₆, λ₄, λ₆ für b = 1.2…3; λ₄-Minimum 4.7799 bei b* = 1.53279; inf λ₄ > 1 vermutet (B6, B7) | numerisch |

Nicht eingetragen: B8 (gp-Fehler beim Einseiten-Test der Tetration an η₃ — technisch).

## Gegenprobe (Pflicht)

| Befund | geprüft wie | Ergebnis |
|---|---|---|
| B5 Sprung an η₅ | `e4c_stencil.py 20 1.737355925151687573101390 6 0.004,0.006,0.008,0.012` — neue Abstände (neue Basen, keine Treffer im Zwischenspeicher), 23 s; Ausgabe `daten/kontrolle_r6_…_0.004,0.006,0.008,0.012_20.json`, `daten/gegenprobe_stencil_out.txt` | R − L = 6.1013·10⁻⁶ bei x = ½ (Bericht 6.10127·10⁻⁶) und 3.6960·10⁻⁶ bei x = 3/2 (Bericht 3.69604·10⁻⁶); mit L′(3/2) = 0.44629 folgt θ(3/2) = 8.2816·10⁻⁶ = θ(½) — Periodizität bestätigt. Einseitige Grenzwerte beider Stencils gleich auf ~10⁻¹⁰, der Sprung ist also kein Fit-Abbruch. **Bestanden.** Die Deutung (echte Kneser-Hexation oder Zwei-Bogen-Auswahl des Motors) bleibt offen — beide Seiten kommen aus demselben Motor. |
| B4 an η₄ | `e4c_stencil.py 20 1.6353244967152763993 5 0.004,0.006,0.008,0.012` — neue Abstände, 2 s | links − rechts 8.7·10⁻¹² (x = ½), 1.5·10⁻¹¹ (x = 3/2) = Fit-Rest bei diesen Abständen; Ableitung in b 0.366123 beidseits (= c₁ = 0.36612345). Sechs Grössenordnungen unter dem Sprung an η₅. **Bestanden.** |
| B6 Landkarte | gegen unabhängige Werte im Graphen | p₄(2) = L₂ = −1.7439091761327, λ₄(2) = 5.1717972635 (kn.beob.pentation_werte); σ₆(1.5) = −3.62649, σ₆(1.7) = −3.69422 (kn.satz.hexation_struktur); p₆(2) = −3.6994078, λ₆(2) = 9.493087 (r18_02_kn); λ₄ bei η₄ und η₅ (4.806, 4.874; kn.beob.tetra_euler, kn.beob.stufe5) liegt glatt zwischen den Rasterwerten. **Bestanden.** λ₄-Minimum nur eine Präzision — so vermerkt. |
| B1 (Beweis) | gelesen | Induktion richtig. Neuheit: **bekannt** (Nixon/Trappmann, Forum tid 1177, Beitrag 2 — schon in der Recherche von kn.satz.hexation_struktur zitiert), nicht „wohl Folklore“. |
| B2 (Beweis) | gelesen, Intervalle nachgerechnet | Richtig unter Monotonie/Stetigkeit. Kleine Berichtigung: „H_r: ℝ → (p_{r−1}, ∞)“ für ungerade r gilt nur für b > η_{r−1}; darunter ist der Wertebereich nach oben beschränkt (durch den kleinsten positiven Fixpunkt von H_{r−1}). Die Aussage links von 0 bleibt unberührt; in der Eintragung so formuliert (nur H_r ↓ p_{r−1}). |
| B3 | gelesen | **Berichtigt:** Für gerade r wechselt H_r die Bauweise an η_{r−1} selbst (Tetration an η₃); die Schranke muss „η_{2j+1} ≤ η_{r−1}“ heissen, nicht „<“. Für ungerade r: η_{2j+1} ≤ η_{r−2}. Als Vermutung eingetragen (Herleitung mit zwei Lücken). |

Kein Befund fiel durch; berichtigt: B1 (Neuheit), B2 (Wertebereich), B3 (Schranke).

## Änderungen am Bestand (Vorschläge)

1. **`kn.konstr.hexation_regulaer`**: `pruefung.f5_neue_frage` 'Geht $H_b$ für $b\uparrow\eta_5$ stetig in eine Kneser-artige Hexation über (wie die Tetration bei η)?' → 'Vorläufig nein: die Motor-Hexation springt an η₅ um eine 1-periodische Höhenverschiebung θ ~ 10⁻⁵ (kn.beob.hexation_sprung_eta5); ob das an der Zwei-Bogen-Auswahl liegt, ist offen.'
2. **`kl.verm.eta_glatt_nicht_analytisch`**: `pruefung.f4_uebertragbar` 'Auf η-Übergänge höherer Ränge (Pentation bei η₅).' → 'An η₅ (Bauwechsel der Hexation) passt das Bild nicht: dort ein Sprung θ ~ 10⁻⁵ (kn.beob.hexation_sprung_eta5); an η₄ ist die Pentation analytisch (kn.beob.pentation_analytisch_eta4).' (Nebenbei: an η₅ wechselt die Hexation, nicht die Pentation.)
3. **`kn.satz.hexation_struktur`**: `pruefung.f4_uebertragbar` (Muster gerade/ungerade) ergänzen: '— für alle Ränge: kn.satz.negative_struktur (bedingt), ganze Werte: kn.satz.ganze_negative_hoehen'.
4. **`kn.konstr.pentation`**: `bedingungen_text` ergänzen: 'λ_b > 1 auf dem Raster b = 1.2…3 gemessen, Minimum 4.7799 bei b* = 1.53279 (kn.beob.landkarte_p4_p6)'.

## Strang und Ideenkarte

- **Neuer Strang** `kn.raenge_basis`, Name „Die Hierarchie als Funktion der Basis“, Mitglieder: die sechs ids oben; `verlauf`: „Links von null wechseln sich die Ränge ab (Pol und Fixpunkt bei geraden, keiner bei ungeraden); nur die ungeraden Grenzbasen η₃, η₅, … sind Bruchstellen — an η₄ ist die Pentation analytisch, an η₅ springt die Hexation um eine periodische Höhenverschiebung θ ~ 10⁻⁵.“
  (Alternative: an `kn.hoehere_stufen` anhängen; der Strang hat aber schon 21 Mitglieder.)
- **Ideenkarte**: neue Unteridee unter `stufen`, id `stufen.basis`, Name „Höhere Ränge als Funktion der Basis“, `kurz`: „Gerade Grenzbasen sind nur Schwellen, ungerade echte Bruchstellen; an η₅ passen reguläre und Kneser-artige Hexation nicht zusammen.“, `straenge: [kn.raenge_basis]`.

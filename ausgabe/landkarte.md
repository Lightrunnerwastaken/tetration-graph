# Landkarte der Tetration

Hauptideen und Unterideen mit Füllstand (Status der Aussagen), Herkunft (Literatur / Forum / nur eigene Arbeit) und den offenen Aussagen (Vermutungen, Fragen, Annahmen, bedingte Sätze; nach Tragweite). Erzeugt von baue.py.

## Fraktionale Fortsetzung (97 Aussagen)

Wie man T(x) zwischen den ganzen Höhen festlegt — und welche Wahl die natürliche ist.

Stand: bewiesen 44, bedingt 9, numerisch 17, Vermutung 14, offen 10, verworfen 3 · Literatur 50 / Forum 17 / eigen 30

- **Halbe Exponentialfunktion** (26): bewiesen 8, bedingt 6, numerisch 6, Vermutung 4, verworfen 2 · Literatur 7 / Forum 1 / eigen 18 · offen: kn.satz.halbexp_beschraenkt, kn.satz.polynom_halbexp_hauptsatz, kn.satz.polynom_halbexp_reduktion, kn.satz.polynom_hauptsatz_kr0, kn.verm.polynom_halbexp, kn.verm.polynom_halbexp_eckmoden, kn.verm.polynom_halbexp_kern, kn.verm.polynom_halbexp_schranke
- **Andere Konstruktionen** (23): bewiesen 10, bedingt 1, numerisch 3, Vermutung 5, offen 4 · Literatur 10 / Forum 11 / eigen 2 · offen: kn.frage.abel_bruchsumme, kn.frage.tpid14, kn.frage.tpid16, kn.frage.tpid17, kn.verm.carleman_eigenwerte, kn.verm.carleman_kneser, kn.verm.szekeres_wachstum, kn.verm.tpid3 …
- **Regulär unter η** (17): bewiesen 9, bedingt 2, numerisch 2, Vermutung 3, offen 1 · Literatur 10 / Forum 3 / eigen 4 · offen: kl.verm.streifenbreite, kl.verm.vollmonoton, kl.frage.tpid5, kl.satz.folgen_vollmonoton, kl.satz.newton_reihe, kl.verm.faktor
- **Der Grenzfall η** (14): bewiesen 8, numerisch 2, Vermutung 1, offen 3 · Literatur 9 / Forum 2 / eigen 3 · offen: br.frage.implosion, kl.frage.eta_analytisch, kl.frage.stetig_eta, kl.verm.lambert_zeta
- **Kneser über η** (9): bewiesen 5, numerisch 3, verworfen 1 · Literatur 7 / Forum 0 / eigen 2
- **Eindeutigkeit** (7): bewiesen 4, numerisch 1, Vermutung 1, offen 1 · Literatur 6 / Forum 0 / eigen 1 · offen: kn.frage.beta_methode, kn.verm.tk_c2
- **Reelle Kennzeichnung?** (1): offen 1 · Literatur 1 / Forum 0 / eigen 0 · offen: br.frage.bohr_mollerup

## Die Gleichung T(x+1) = aᵀ⁽ˣ⁾ (70 Aussagen)

Was jede Lösung der Grundgleichung erfüllt, ganz gleich, wie man sie zwischen den ganzen Höhen festlegt.

Stand: bewiesen 58, numerisch 7, Vermutung 4, offen 1 · Literatur 20 / Forum 19 / eigen 31

- **Werte, Ableitungen, Gestalt** (28): bewiesen 20, numerisch 5, Vermutung 3 · Literatur 4 / Forum 10 / eigen 14 · offen: gr.verm.ein_wendepunkt, gr.verm.tpid7, gr.verm.ungerade_ableitungen
- **Integrale** (13): bewiesen 9, numerisch 2, Vermutung 1, offen 1 · Literatur 2 / Forum 1 / eigen 10 · offen: gr.frage.integral_form, gr.verm.turmintegral_rate
- **Begriffe** (8): bewiesen 8 · Literatur 6 / Forum 0 / eigen 2
- **Freiheit der Lösungen** (6): bewiesen 6 · Literatur 3 / Forum 1 / eigen 2
- **Was fest ist, was frei ist** (5): bewiesen 5 · Literatur 2 / Forum 2 / eigen 1
- **Superlogarithmus** (5): bewiesen 5 · Literatur 1 / Forum 2 / eigen 2
- **Superwurzeln** (5): bewiesen 5 · Literatur 2 / Forum 3 / eigen 0

## Basiswechsel T ∘ slog (69 Aussagen)

Was entsteht, wenn man die Tetration zu einer Basis mit dem Superlogarithmus einer anderen zusammensetzt.

Stand: bewiesen 19, bedingt 17, numerisch 11, Vermutung 6, offen 6, verworfen 10 · Literatur 7 / Forum 6 / eigen 56

- **Antwortkern und Basisableitung** (18): bewiesen 3, bedingt 5, numerisch 3, Vermutung 2, verworfen 5 · Literatur 1 / Forum 2 / eigen 15 · offen: bw.ann.K1K2, bw.satz.kern_eindeutig, bw.ann.basis_diffbar, bw.satz.kern_existenz, bw.prop.hub, bw.satz.antwortkern
- **Moden und Singularitäten** (13): bewiesen 2, bedingt 1, numerisch 5, Vermutung 2, offen 1, verworfen 2 · Literatur 1 / Forum 2 / eigen 10 · offen: bw.frage.streifen_gevrey, bw.verm.nirgends_analytisch, bw.verm.streifen
- **Phasengesetz und Welle** (10): bewiesen 2, bedingt 3, numerisch 2, Vermutung 1, offen 2 · Literatur 3 / Forum 1 / eigen 6 · offen: bw.satz.phasengesetz, bw.frage.halbiter_vergleich, bw.frage.nixon_eta_basiswechsel, bw.lem.offset, bw.verm.adressen
- **Gleicher Fixpunkt** (8): bewiesen 4, bedingt 1, offen 2, verworfen 1 · Literatur 2 / Forum 0 / eigen 6 · offen: kl.frage.basiswechsel_kneser, kl.frage.gleiches_L_kneser
- **Log-Peeling und Grenzkanal** (8): bewiesen 6, bedingt 1, Vermutung 1 · Literatur 0 / Forum 1 / eigen 7 · offen: bw.ann.hyp22, bw.satz.schwanzkonvergenz
- **Kreisabbildung** (6): bewiesen 1, bedingt 4, verworfen 1 · Literatur 0 / Forum 0 / eigen 6 · offen: bw.prop.kozykel, bw.satz.homoeo, bw.satz.antisymmetrie, bw.satz.c1_produkt
- **Rekonstruktion und Zertifikat** (6): bewiesen 1, bedingt 2, numerisch 1, offen 1, verworfen 1 · Literatur 0 / Forum 0 / eigen 6 · offen: bw.frage.intervallmarsch, bw.satz.anker_rekonstruktion, bw.zert.8r1

## Höhere Stufen (49 Aussagen)

Pentation, Hexation und die ganze Leiter der Hyperoperationen.

Stand: bewiesen 13, bedingt 7, numerisch 19, Vermutung 4, offen 4, verworfen 2 · Literatur 11 / Forum 16 / eigen 22

- **Pentation und Hexation** (17): bewiesen 1, bedingt 7, numerisch 7, Vermutung 1, offen 1 · Literatur 4 / Forum 7 / eigen 6 · offen: kn.konstr.pentation, kn.satz.hexation_eindeutig_g0, kn.frage.hexation_kneser, kn.konstr.hexation_regulaer, kn.satz.hexation_gebiete, kn.satz.hexation_keine_halbebene, kn.satz.hexation_logpunkt, kn.satz.hexation_struktur …
- **ηₙ → 2 und eₙ → 4** (17): bewiesen 5, numerisch 7, Vermutung 3, offen 2 · Literatur 1 / Forum 9 / eigen 7 · offen: kn.frage.andere_etas, kn.frage.pentation_untergrenze, kn.verm.beschraenkt_unter_2, kn.verm.tpid10, kn.verm.tpid11
- **Gebrochene Ränge** (11): bewiesen 3, numerisch 5, offen 1, verworfen 2 · Literatur 2 / Forum 0 / eigen 9 · offen: kn.frage.halber_rang
- **Hyperoperationen** (4): bewiesen 4 · Literatur 4 / Forum 0 / eigen 0

## Unendlicher Turm und Basen-Ebene (43 Aussagen)

Für welche Basen der unendlich hohe Turm einen Grenzwert hat, und welche Bilder dabei in der komplexen Ebene entstehen.

Stand: bewiesen 30, bedingt 1, numerisch 6, Vermutung 1, offen 3, verworfen 2 · Literatur 27 / Forum 5 / eigen 11

- **Regionen und Fraktale** (14): bewiesen 9, numerisch 4, Vermutung 1 · Literatur 9 / Forum 0 / eigen 5 · offen: fp.verm.knospen_drehzahl
- **Konvergenz des Turms** (13): bewiesen 11, numerisch 1, offen 1 · Literatur 9 / Forum 4 / eigen 0 · offen: fp.frage.tpid12_13
- **Spiralen und Mittelwerte** (10): bewiesen 4, bedingt 1, numerisch 1, offen 2, verworfen 2 · Literatur 4 / Forum 0 / eigen 6 · offen: fp.frage.holomorphe_familie, fp.verm.doppelmittel
- **Periode 2 unter e⁻ᵉ** (6): bewiesen 6 · Literatur 5 / Forum 1 / eigen 0

## Zahlen und Arithmetik (30 Aussagen)

Welche Werte rational, algebraisch oder transzendent sind, und Türme modulo m.

Stand: bewiesen 20, bedingt 3, numerisch 1, Vermutung 3, offen 3 · Literatur 19 / Forum 2 / eigen 9

- **Transzendenz** (23): bewiesen 15, bedingt 3, Vermutung 2, offen 3 · Literatur 15 / Forum 2 / eigen 6 · offen: ar.verm.schanuel, ar.frage.e_hoch_e, ar.frage.ef_fortsetzung, ar.frage.ganze_tuerme, ar.satz.schanuel_orbit, ar.verm.gelfond_turm
- **Kongruenzen** (4): bewiesen 4 · Literatur 4 / Forum 0 / eigen 0
- **Kneser–Schanuel** (3): bewiesen 1, numerisch 1, Vermutung 1 · Literatur 0 / Forum 0 / eigen 3 · offen: ar.verm.ks

## Dynamik von bᶻ und der Fixpunkt L (24 Aussagen)

Was passiert, wenn man b^z immer wieder anwendet — und welche Rolle die Fixpunkte spielen.

Stand: bewiesen 24 · Literatur 21 / Forum 1 / eigen 2

- **Parabolisch** (10): bewiesen 10 · Literatur 10 / Forum 0 / eigen 0
- **Anziehend und abstossend** (8): bewiesen 8 · Literatur 8 / Forum 0 / eigen 0
- **Fixpunkte von bᶻ** (5): bewiesen 5 · Literatur 2 / Forum 1 / eigen 2
- **Indifferent** (1): bewiesen 1 · Literatur 1 / Forum 0 / eigen 0

## Brücken (11 Aussagen)

Wo die Tetration in Logik, Wachstumstheorie und Fourier-Analysis auftaucht.

Stand: bewiesen 8, offen 2, verworfen 1 · Literatur 10 / Forum 0 / eigen 1

- **Logik und Wachstum** (7): bewiesen 5, offen 2 · Literatur 7 / Forum 0 / eigen 0 · offen: br.frage.kneser_hardy, br.frage.transexp
- **Kreis und Fourier** (4): bewiesen 3, verworfen 1 · Literatur 3 / Forum 0 / eigen 1

## Stärkste Kopplungen zwischen Hauptideen

- Fraktionale Fortsetzung — Die Gleichung T(x+1) = aᵀ⁽ˣ⁾: 52.04
- Basiswechsel T ∘ slog — Fraktionale Fortsetzung: 40.05
- Dynamik von bᶻ und der Fixpunkt L — Fraktionale Fortsetzung: 29.06
- Fraktionale Fortsetzung — Höhere Stufen: 28.81
- Die Gleichung T(x+1) = aᵀ⁽ˣ⁾ — Unendlicher Turm und Basen-Ebene: 22.61
- Basiswechsel T ∘ slog — Die Gleichung T(x+1) = aᵀ⁽ˣ⁾: 21.35
- Die Gleichung T(x+1) = aᵀ⁽ˣ⁾ — Höhere Stufen: 18.37
- Fraktionale Fortsetzung — Unendlicher Turm und Basen-Ebene: 17.25
- Dynamik von bᶻ und der Fixpunkt L — Unendlicher Turm und Basen-Ebene: 14.22
- Basiswechsel T ∘ slog — Unendlicher Turm und Basen-Ebene: 12.32

## Stärkste Kopplungen über Ideengrenzen

- Anziehend und abstossend (Dynamik von bᶻ und der Fixpunkt L) — Regulär unter η (Fraktionale Fortsetzung): 8.14
- Parabolisch (Dynamik von bᶻ und der Fixpunkt L) — Der Grenzfall η (Fraktionale Fortsetzung): 7.77
- Gleicher Fixpunkt (Basiswechsel T ∘ slog) — Anziehend und abstossend (Dynamik von bᶻ und der Fixpunkt L): 5.67
- Andere Konstruktionen (Fraktionale Fortsetzung) — Werte, Ableitungen, Gestalt (Die Gleichung T(x+1) = aᵀ⁽ˣ⁾): 5.61
- Regulär unter η (Fraktionale Fortsetzung) — Werte, Ableitungen, Gestalt (Die Gleichung T(x+1) = aᵀ⁽ˣ⁾): 5.61
- Werte, Ableitungen, Gestalt (Die Gleichung T(x+1) = aᵀ⁽ˣ⁾) — Konvergenz des Turms (Unendlicher Turm und Basen-Ebene): 4.6
- Fixpunkte von bᶻ (Dynamik von bᶻ und der Fixpunkt L) — Regionen und Fraktale (Unendlicher Turm und Basen-Ebene): 4.5
- Antwortkern und Basisableitung (Basiswechsel T ∘ slog) — Freiheit der Lösungen (Die Gleichung T(x+1) = aᵀ⁽ˣ⁾): 4.4
- Gleicher Fixpunkt (Basiswechsel T ∘ slog) — Regulär unter η (Fraktionale Fortsetzung): 4.14
- Kneser über η (Fraktionale Fortsetzung) — Pentation und Hexation (Höhere Stufen): 4.04
- Andere Konstruktionen (Fraktionale Fortsetzung) — Pentation und Hexation (Höhere Stufen): 3.97
- Antwortkern und Basisableitung (Basiswechsel T ∘ slog) — Werte, Ableitungen, Gestalt (Die Gleichung T(x+1) = aᵀ⁽ˣ⁾): 3.87
- Regulär unter η (Fraktionale Fortsetzung) — Begriffe (Die Gleichung T(x+1) = aᵀ⁽ˣ⁾): 3.6
- Anziehend und abstossend (Dynamik von bᶻ und der Fixpunkt L) — Begriffe (Die Gleichung T(x+1) = aᵀ⁽ˣ⁾): 3.43
- Eindeutigkeit (Fraktionale Fortsetzung) — Freiheit der Lösungen (Die Gleichung T(x+1) = aᵀ⁽ˣ⁾): 3.42

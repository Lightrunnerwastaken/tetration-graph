# Prototypen: Kneser-artige Hexation über η₅ (09.10.2026)

Nicht als Experimente registriert — Arbeitsstand zu `kn.frage.hexation_kneser`. Die Skripte erwarten ihre Nachbarn im selben Ordner
(Pfade in `sys.path.insert` auf den alten Arbeitsordner zeigen; vor einer Wiederaufnahme anpassen).

| Datei | Inhalt | Ergebnis |
|---|---|---|
| `hexa.py` | reguläre Superfunktion am **anziehenden** Fixpunkt p der Pentation (Basis 2) über P⁻¹ | springt — erklärt durch `kn.satz.kein_kneser_anziehend` |
| `hexq.py` | reguläre Superfunktion F_q am **abstossenden** Fixpunkt q = −1.99295 + 1.44513i (nur Vorwärtsschritte) | glatt |
| `hex_lauf.py` (+ `.json`, `.log`) | Kurve C (F_q reell) als Graph g(t), Theodorsen, Normierung | 2↑↑↑↑1.5 ≈ 2.6760 / 2.6821 / 2.6736 (K = 32/64/128), Theodorsen ab K = 64 nicht konvergent |
| `hex_fornberg.py` | Kurve über den reellen Wert w, Fornberg | Kurve gut, Gauss–Newton divergiert |
| `hex_kurve.py`, `hex_ecke.py` | |F_q′| entlang C; Newton auf F_q′ | kritischer Punkt ζ_c ≈ 0.4031293 − 1.7097883i, Wert ≈ −0.7715464141 (dicht bei P(L_T) = −0.7715464478), Konvergenz nur linear |
| `hex_kouz.py`, `hex_zweig.py` | Kouznetsovs Cauchy-Methode im Streifen |Re z| ≤ 1 | divergiert; Startwert H(iy) springt bei y ≈ 1.5 (Hauptzweig der Kneser-Tetration in der Kette gekreuzt) |

| `blatt.py` | Blattverfolgung: T und P entlang von Wegen fortgesetzt (Zweigzahlen nach Stetigkeit, adaptive Verfeinerung `verfolge`); riesige Zwischenwerte als Abstand zum Hauptzweig mit exaktem Logarithmus ln b·T(x − 1) | `test_blatt.py`: reell = Hauptzweig (5e−27), Kreis kehrt zurück, Weg durch Werte ~5000 = Hauptzweig; `test_blatt_kontrolle.py`: P(1 + 3i) verfolgt = 1.634077394 + 0.960120499i auf drei Arten gleich |
| `fixpunkt_verfolgt.py` | P entlang gerader und „oben herum“-Wege bis q | q bleibt Fixpunkt der fortgesetzten Pentation, μ unverändert |
| `start_asymptotisch.py`, `newton_hex.py` | Startwert (grob bis y = 1.5, darüber q + C·e^{−1.068(y − 1.5)}); Newton für H(iy) im Rechteck [0, 2] × [−A, A] mit Cauchy-Gleichung auf Re z = 1 | aufgegeben: P(P(H₀)) läuft bei y ≈ 2.845 fast genau durch einen Verzweigungspunkt; grundsätzlich ungültig, weil H im Rechteck verzweigt ist |
| `diag_weg.py`, `diag_sprung.py`, `diag_seite.py` | wo der Weg nicht auflösbar ist; welche Stufe der Kette; Umweg links/rechts | Stufe 5 der Kette trifft −2; Hauptzweig-P hat einen Schnitt auf Im w = π/ln λ_T ab w_b ≈ −1.979 + 1.912i; beide Umwege enden nicht in q |
| `karte_P.py` (+ `.png`) | Phasenbild der Hauptzweig-Pentation | periodisch mit 2πi/ln λ_T = 3.824i, Schnitt ab w_b nach rechts |
| `karte_H.py` (+ `_64.png`) | Phasenbild von H aus dem groben θ (Hauptzweig) | oberhalb Höhe 3 glatt, Nullstelle-Kopie bei ≈ 1.95 + 2.97i; darunter Hauptzweig-Sprünge (nicht verwertbar) |
| `monodromie.py`, `exponent.py` | Umläufe mit Blattverfolgung um ζ_c − 2 + m + Π | m = 0: +9.06472i = 2πi/ln 2; m = 1, 2: ≠ 0; Kontrollen: 0; Exponent α ≡ 0.37070 (vorhergesagt ln λ_T/ln 2 = 2.370666) → `kn.beob.hexation_verzweigt` |
| `cache_drossel.py` | gemeinsamer Zwischenspeicher nur lesen, neue Werte nach `kneser_cache_proto.json` (alle 5 min und beim Beenden) | ohne das schrieb jeder Motoraufruf ~190 MB |

| `probe_fork.py`, `probe_radius.py`, `probe_b.py` | Kneser-Hexation aus dem gp-tetration-Fork (knhexinit, exp-079) | Basis 2 in 12 s mit 32.6 Stellen, 2↑↑↑↑1.5 = 2.672893959342619489…; θ-Konstante −c₀ = 5.16210 + 1.42703i = unser XS + θ̂₀ + Π bis 0.003 → Verzweigungspunkte genau bei −0.74076 + m + 2.74581i; Abel-Reihe um Re q hat ihre nächste Singularität bei q (1.01 Im p), nicht bei „0.64 Im p“ |
| `eta6.py` | η₆, e₆ über knhexinit: Gitter + Newton auf H′ = 1, Sekanten in b | dps 30 und 45 gleich auf 32 Stellen: η₆ = 1.79912123612881268672356700702897599, e₆ = 3.44223464771689323 → `e39_hexation_kneser.py` |

| `probe_egrenze.py`, `probe_reste.py`, `probe_ecke.py`, `probe_pentgenau.py`, `probe_luecke.py`, `probe_fernbilder.py` | e-Grenze des Fork-Verfahrens (~26 Stellen bei e) | Ecken liegen in W unter hmin (e 1.412, 2 1.231) → θ-Abbruch harmlos; pent/ipent auf 1e-48…1e-55 genau; Lücke φ < 0.05 harmlos; Ursache: FE-Gleichungen mit Bild bis 1.22 Bogenradien draussen — ohne die jenseits 1.05: 25.7 → 34.1 Stellen (dps 45; bei ≤ 1.0 kippt es auf 18). Mit nach rechts verschobenem knc (exp-081-Zwischenstand) begrenzt eine Singularität bei 0.33 ± 1.37i (Abstand 1.46 statt |c − q| = 2.5) |

| `karte_reell.py`, `karte_unten.py` (+ `.png`, `.npz`) | reelle Menge von F_q mit Blattverfolgung: Spalten von oben bzw. Zeilen von links | über C bis C + Π und unter C bis C − Π keine weitere verschiebungsfeste reelle Kurve; die Kneser-Kurve C = Bögen von Ecke zu Ecke ab dem Log-Punkt ζ_c − 2; Schleifen unter den Ecken nur auf dem Blatt „von oben“ → `kn.verm.hexation_streifen` |
| `streifen.py` (+ `.png`) | H auf den Rändern von S = [−1.8, −0.8] (θ grob, K = 128, F_q verfolgt) | Randbilder einfach und disjunkt → H auf S injektiv; w_b im Bild; Streifen ohne Singularität: s_b < Re z < −0.741 (Breite 3.03 ≈ Re Π) |

| `gebiet_g0.py` (+ `.png`, `.log`) | festes Fundamentalgebiet G₀ zwischen γ₀ = [q̄, q] (senkrecht durch Re q) und P(γ₀); α_K auf γ₀ unten über ihex, oben über θ-Darstellung | Überlapp 1e-14…1e-18; γ₀ ∩ P(γ₀) = {q, q̄}; max Re α_K(γ₀) = −2.994 → H_K bildet den Streifen biholomorph auf G₀ ab → `kn.satz.hexation_eindeutig_g0` |

| `gebiete_gc.py` (+ `gebiete_gc*.log`) | Kreisbogen-Familie G_c (γ_c von q̄ über c nach q): Jordan-Rand? Bildstreifen zwischen α_K(γ_c) und +1 frei von Verzweigungspunkten z_{m,k}? | c ∈ {−2.5, Re q, −1.6, −1.2, −0.9, −0.6, −0.3, 0, 0.4, 1}: ja → Eindeutigkeit auf allen (`kn.satz.hexation_gebiete`); c = 2, 3: P(γ_c) schneidet γ_c (kein Fundamentalgebiet). ihex mit `alarm(2, …)` begrenzen (gibt bei Zeitablauf t_ERROR zurück, nicht abfangbar mit iferr) |

| `umlauf_p.py`, `fortsetzung_p.py` | α_K auf Kreisen um den anziehenden Fixpunkt p (ihex je Punkt); Newton-Fortsetzung von α_K über hex | ihex deckt die Kreise nur teilweise; hex reicht nur bis Im z ≈ 0.7, Fortsetzung bricht bei 0.61 + 0.79i ab |
| `umlauf_w.py`, `umlauf_w2.py` | Umlauf der q-Abel-Koordinate W = α_q (knwalk bzw. knalphap je Punkt) um q, um einen gewöhnlichen Punkt und um p | um q: ΔW = 1.000000·Π; gewöhnlich: 0; um p: Sprünge ≈ ±Π_p (+ Bruchteile von Π), nicht durch Vielfache von Π behebbar → α_q ist an p verzweigt |
| `karte_rechts.py` (+ `.png`, `.npz`) | F_q als Hauptkette des Forks (pent^n(q + Ψ(μ^(ζ−n))), Ψ zur Laufzeit aus ksig) bis Re ζ = 7 | im Bereich von karte_reell gleich der Blattverfolgung auf 1e−14, 1 min statt 12; zwischen C und Im 0.2 bis Re 4.7 nur C; oberhalb (rechts von 1.41) Schnitt der Hauptkette |
| `logpunkte.py`, `strahl.py` (+ `.log`, `strahl_*.npy`) | Log-Punkte Z = ζ* + 1 mit F_q(ζ*) = w_b + jΠ_T (Newton-Gitter, j = −1 … 2) im Band zwischen C und C + Π; reeller Strahl per Newton-Fortsetzung bis L_T | vier Log-Punkte: ζ_c − 2, ζ_c − 2 + Π, 2.68744 + 0.99506i, 3.63707 + 0.95142i; nur der Kneser-Strahl schliesst bei Z + 1 → `kn.satz.hexation_logpunkt`; ζ_c = 0.40306781 − 1.70988020i |

Eindeutigkeit (09.10.2026): Halbebenen-Kennzeichnung à la Paulsen–Cowgill leer (`kn.satz.hexation_keine_halbebene`, bedingt auf die Log-Verzweigung);
mit festem Fundamentalgebiet G₀ eindeutig (`kn.satz.hexation_eindeutig_g0`: Eindeutigkeit bewiesen, Existenz numerisch); offen nur die
Gebietsunabhängigkeit (`kn.verm.hexation_streifen`). Hinweis: Im Fork exp-081 bricht `knforce=[0.6,0,0,0,0]` (alter Bogen um Re q) ab, weil die
neue Bildpunkt-Regel ihn verwirft; `ihex` hängt für Punkte von γ₀ oberhalb ~0.3 Im q — dort die θ-Darstellung nehmen.

Die Probe-Skripte lesen den Fork nur (FatouGP.eval_batch; kn2solve ist eine Laufzeit-Kopie von knsolve). Der Fork wurde am 09.10. ab ~14:40 von der Sitzung „gp-tetration Optimierungsstand“ umgebaut (exp-081, Mittelpunkt und Radius frei) — Ergebnisse hängen vom Fork-Stand ab.

Stand 09.10.2026: Das θ-Verfahren mit Eckenbasis ist nicht mehr nötig — das Fork-Verfahren rechnet die Abel-Funktion nahe q und an der Achse über die
Funktionalgleichung und umgeht die Ecken. Offen: Eindeutigkeit (H ist in der oberen Halbebene verzweigt), Übergang bei η₅, die e-Grenze (~26 Stellen) des Fork-Verfahrens.
Vergleichswert: S. Levenstein (Tetration Forum tid 1263, 2020): 2↑↑↑↑1.5 ≈ 2.6729.

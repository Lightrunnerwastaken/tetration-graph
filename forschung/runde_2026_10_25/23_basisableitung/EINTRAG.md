# Eintrag Bericht 23 (Basisableitung: Vorzeichen und Formeln) in den Graphen

Stand 25.10.2026 17:55. Dateien: `daten/aussagen/r25_23_gr.yaml`, `daten/aussagen/r25_23_kl.yaml`; Theorie: `daten/theorie/r25_23_gr.yaml`, `daten/theorie/r25_23_kl.yaml`.
`python baue.py --nur-pruefen`: 0 Fehler in diesen Dateien (nur die erwarteten Hinweise).

## Angelegte ids
- gr.satz.elastizitaet_stufe | E(x+1) = (1 + E(x)) ln T(x+1), E(0) = 0, E(1) = ln b, Summenform an ganzen Höhen (B1) — bewiesen
- gr.satz.vorzeichen_reduktion | Vorzeichenmuster ⇔ max(−1, ln T) < E < 0 auf (−1, 0); slog-Form; Schnittpunkte; D → +∞ bei −2 (B2, B7, B16) — bewiesen
- gr.satz.elastizitaet_rand | E(−1⁺) = ε₀ − 1, E′(0) = ε₀T′(0); Muster verlangt 0 ≤ ε₀ ≤ 1 (B3) — bewiesen
- gr.beob.basis_monoton | Muster an 13 Basen (fein) und 39 Basen (dicht), E(−½) über die Basen (B4–B6, B8) — numerisch
- gr.beob.eps0_tabelle | Tabelle ε₀(b), T′(0) von 1 + 10⁻⁶ bis 10⁴; Steigung für grosse b (B9, B11) — numerisch
- gr.beob.wackeln_basen | u = u₀ + w; Spanne von w über die Basen, 0.5–0.64 % (B12–B14) — numerisch (Struktur bewiesen)
- gr.frage.vorzeichen_kneser | Beweis der inneren Ungleichungen; ε₀ für b → ∞ (N2, Offene Fragen 1, 2, 5) — offen
- kl.satz.eps0_b_gegen_1 | ε₀ = 1 − 1/ln(1/ln b) + Rest, T′(0) ≈ (a − 1) ln(1/λ) (B10) — skizze
- kl.satz.kanalgesetz_vorwaertssumme | S_b(0) ~ (π/2)√(2η/e)(b − η)^{−3/2}, Korrektur 1 − 0.644 δ (B15) — skizze

Nicht eigens eingetragen: B7 und B16 (in gr.satz.vorzeichen_reduktion), B11 (in gr.beob.eps0_tabelle), B14 (in gr.beob.wackeln_basen), N1 (Motorgrenze, in den Prüffragen von gr.beob.eps0_tabelle).
Gebiet geändert gegenüber dem Vorschlag im Bericht: kn.beob.basis_monoton → gr.beob.basis_monoton, kn.beob.eps0_tabelle → gr.beob.eps0_tabelle, gr.satz.kanalgesetz_… → kl.satz.kanalgesetz_vorwaertssumme (damit sie in die Stränge gr.basis bzw. kl.eta passen).

## Gegenprobe
| Befund | geprüft wie | Ergebnis |
|---|---|---|
| B1 | Rechnung nachvollzogen; gleichwertig zu gr.satz.basis_gleichung geteilt durch T(x+1)/(b ln b); Summenform für n = 1, 2 nachgerechnet | bestanden → bewiesen |
| B2 (+B7, B16) | Fallunterscheidung je Streifen nachgeprüft (Vorzeichen von T, ln T); slog-Form aus D + T′∂_b slog = 0; Rand bei −2 aus der umgestellten Stufengleichung | bestanden → bewiesen; Voraussetzung T′ > 0 für die slog-Form ergänzt |
| B3 | l'Hôpital bei −1, ∂_b ln ln b = 1/(b ln b), E′(0) nachgerechnet; Notwendigkeit 0 ≤ ε₀ ≤ 1 aus erster Ordnung | bestanden → bewiesen |
| B4, B8 | neue Rechnung an einer **neuen Basis**: `d2_elastizitaet.py 28 9 9 fein 3 13/10`, dann `a2_muster.py` | b = 3 und 1.3: Muster auf dem ganzen Gitter, E steigend, b = 3: min(E+1) = 0.559, min((E − ln T)/\|ln T\|) = 0.477, E(−½) = −0.2604 (liegt zwischen 2 und 10) |
| B9 | neue Einstellung `d2_elastizitaet.py 32 10 10 rand 2 "exp(exp(-1))+1/1000"` (frische Motorwerte) | ε₀(2) = 0.571578377725955, ε₀(η+10⁻³) = 0.612836763472541, ε₀(1.3) = 0.63440965531712 — auf alle angegebenen Stellen gleich |
| B10 | Herleitung nachvollzogen (Linearisierung, Ableitung nach ln ln b); Formelwerte bei 1 + 10⁻⁶ und 1 + 10⁻⁴ nachgerechnet | Leitterm stimmt; Rest wandert → skizze, Fehlerordnung als offen markiert |
| B12, B13 | Teleskopsumme geprüft; b = 3 aus der neuen Rechnung | Spanne 2.34·10⁻⁴ = 0.58 % von \|u(−½)\|, Periodizität 7.9·10⁻¹⁷ — passt zur Reihe |
| B15 | Herleitung nachgerechnet (α = e²δ/η, β = 1/(2e), ∫dy/q² = π/(2α^{3/2}β^{1/2}), Vorfaktor e²/η ⇒ 1.61946527); `a3_kanal.py` auf die frische Rechnung bei η + 10⁻³ | Verhältnis 0.9993567761 (Bericht 0.99936) — bestanden; Lücke (Kanalnäherung) bleibt → skizze. Nebenbei exakt: S_b(0) = S^Kobi_b(1)/T′_b(0) |

## Änderungen am Bestand
1. `gr.satz.basis_gleichung` (00_grundlagen.yaml), Feld `formen`: neue Form anhängen
   `{tex: 'E(x+1)=(1+E(x))\ln T(x+1),\ E=b\ln b\,D/T', notiz: 'logarithmische Form (gr.satz.elastizitaet_stufe)'}`. Begründung: dieselbe Gleichung in der Form, die das Vorzeichen durchschaubar macht.
2. `gr.beob.kneser_pde_wackeln` (00_grundlagen.yaml), Feld `pruefung.f4_uebertragbar`: alt „Jede Basis über η; die Spanne fällt von 2 zu e.“ → neu
   „Jede Basis über η; die Spanne fällt von 1.5·10⁻³ (knapp über η) über 5.7·10⁻⁴ (2) bis 7·10⁻⁹ (10⁴), relativ zu u von 0.64 % auf 0.50 % (gr.beob.wackeln_basen).“ Begründung: jetzt über sieben Basen gemessen.
3. `kn.beob.grosse_basis` (r18_03_kn.yaml), Feld `pruefung.f5_neue_frage`: alt „Strebt x_w(b) gegen eine Grenze (−3/2? −2?) und wächst T′(0) wie ln ln b?“ → neu
   „Strebt x_w(b) gegen eine Grenze (−3/2? −2?)? T′(0) wächst zwischen 10 und 10⁴ etwa wie (ln b)^{0.45} (ε₀ = ∂ln T′(0)/∂ln ln b ≈ 0.41–0.49, gr.beob.eps0_tabelle), nicht wie ln ln b.“ Begründung: Teilantwort aus der ε₀-Tabelle.

## Vorschlag Stränge
- Bestehender Strang `gr.basis` („Ableitung nach der Basis“): aufnehmen gr.satz.elastizitaet_stufe, gr.satz.vorzeichen_reduktion, gr.satz.elastizitaet_rand, gr.beob.basis_monoton, gr.beob.eps0_tabelle, gr.beob.wackeln_basen, gr.frage.vorzeichen_kneser.
  Neue `verlauf`-Zeile: „Bei ganzer Höhe eine einfache Rekursion; am unendlichen Turm eine Quadratwurzel bei η. Die Elastizität gehorcht E(x+1) = (1 + E) ln T(x+1); damit hängt das gemessene Muster — eine grössere Basis hebt die Kurve, ausser zwischen −1 und 0 — nur an zwei Ungleichungen auf (−1, 0) und an der Randzahl ε₀ ∈ (0, 1).“
- Bestehender Strang `kl.eta`: aufnehmen kl.satz.kanalgesetz_vorwaertssumme; `verlauf` ergänzen um „… knapp über η wächst die Vorwärtssumme der Basisableitung wie (b − η)^{−3/2}.“
- Bestehender Strang `kl.regulaer`: aufnehmen kl.satz.eps0_b_gegen_1 (verlauf unverändert).

## Vorschlag Ideenkarte
- Unteridee `gleichung.werte` („Werte, Ableitungen, Gestalt“, trägt gr.basis): `kurz` ergänzen um „… und wie eine grössere Basis die Kurve hebt und senkt“.
- kl-Einträge laufen über `fortsetzung.eta` und `fortsetzung.regulaer` mit; keine neue Idee nötig.

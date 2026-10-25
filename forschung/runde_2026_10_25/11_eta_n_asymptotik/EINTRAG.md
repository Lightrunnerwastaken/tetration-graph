# Eintrag Bericht 11 (Asymptotik von η_r, e_r; β_r, ι_r) — 25.10.2026

Dateien: `daten/aussagen/r25_11_kn.yaml`, `daten/theorie/r25_11_kn.yaml`. `python baue.py --nur-pruefen`: 0 Fehler, keine Warnung zu diesen ids.

## Angelegte ids

| id | Titel kurz | Status |
|---|---|---|
| kn.satz.kritische_rate_log | e_r − H_{r+1}(η_r,k) = 1/(β_r(k + ι_r ln k + C_r + o(1))) (B3) | bewiesen |
| kn.beob.iota_raenge | ι₃ = 1/3, ι₄…ι₈ ≈ −0.39 … −0.61, Minimum bei r = 7; schrumpfende Radien um e_r (B1, B2, B12) | numerisch |
| kn.beob.beta_raenge | β_r für r = 3…8, wachsend mit fallenden Zuwächsen, Grenzwert offen (B4, B5) | numerisch |
| kn.verm.eta_r_log | 2 − η_r ~ C/(r ln r), C ≈ 2.20; Exponent 1.46 nur Übergang; Test η₉ ≈ 1.88883 gegen 1.88914 (B6–B9) | vermutung |
| kn.gescheitert.e_potenzgesetz | e_r folgt keinem einfachen Gesetz, freie Grenzfits über 4 (B10, B11) | gescheitert |

Nicht eingetragen: B13 (η₉ im Zeitbudget nicht rechenbar — nur als Hinweis in `kn.verm.eta_r_log`), B14 (technischer Abbruch der Oktation bei Radius 0.1).

## Gegenprobe (Pflicht)

| Befund | geprüft wie | Ergebnis |
|---|---|---|
| B1/B4 (ι_r, β_r) | `jets.py --rank 6 --dps 28 --M 24 --rho 0.05` — dritte, neue Cauchy-Einstellung (10.6 s), Ausgabe `cache/jet_r6_d28_M24_rho0.05.json` | β₆ = 0.4429294090481398465, ι₆ = −0.6006192688985848411 — gleich den beiden Berichtseinstellungen auf alle 19 Stellen. **Bestanden.** |
| B4 (β_r) | gegen unabhängige Werte im Graphen | β₄ = 0.29348 (kn.beob.tetra_euler), β₆ = 0.44293 (kn.beob.stufe6), β₇ = 0.4961476452, β₈ = 0.5368242374 (kn.satz.kritische_rate): gleich. Neu ist nur β₅ = 0.37698565404. |
| B6/B7/B9 (v_r, Aitken, Exponenten) | `analyse.py` erneut; dazu eigene Rechnung aus den η-Werten des Graphen (mpmath) | v_r = 1.8303, 2.0222, 2.1135, 2.1596, 2.1827, 2.1935; Aitken 2.1965, 2.2062, 2.2061, 2.2028; lokale Exponenten 1.4619, 1.4708, 1.4705, 1.4662, 1.4603 — wie im Bericht. **Bestanden.** Bemerkung: (2 − η_r)·r^1.5 wächst nahezu linear (2.886 … 2.984), daher die instabilen Aitken-Werte dort; die Zuwächse von v_r halbieren sich fast (Quotienten 0.48, 0.50, 0.50, 0.46). Für C/(r ln r) wäre der lokale Exponent ≈ 1 + 1/ln r (1.48 bei r = 8) — mit den Messwerten verträglich. |
| B3 (Beweis) | Beweis gelesen, Rechnung nachvollzogen | Richtig: −1/g = −1/w + 1 + (c−1)w + O(w²), log(−g) = log(−w) + w + O(w²), also Φ(g) − Φ(w) = 1 + O(w²); ι = 1 − a₃/a₂² ist Écalles résidu itératif; r = 3 gibt ι = 1/3 (a₂ = 1/(2e), a₃ = 1/(6e²)). **Präzisiert:** gebraucht wird die Taylor-Entwicklung bis O(w⁴), also C⁴ (nicht nur C³ wie in kn.satz.kritische_rate) — so in den Voraussetzungen vermerkt. Eingetragen als `bewiesen`, Neuheit `folklore` (klassische Methode). |
| B11 | gegen den Bestand | (4 − e_r)/(2 − η_r) = 2.308 … 2.877 steht schon in `kn.gescheitert.konvexitaetsweg` — nicht neu; nur als Verweis in `kn.gescheitert.e_potenzgesetz`. |
| B12 | gelesen | Nur effektive Radien (k ≤ 8, bei r = 8 k ≤ 5); als Nebensatz in `kn.beob.iota_raenge`, ohne Asymptotik-Anspruch. |

Kein Befund fiel durch. Genauigkeit von ι₈/β₈: begrenzt durch η₈, e₈ (13 bzw. 11 Stellen), nicht durch die Cauchy-Abtastung — so vermerkt.

## Änderungen am Bestand (Vorschläge)

1. **`kn.beob.eta7_eta8`**: `pruefung.f5_neue_frage` 'Warum sinkt der lokale Exponent? Ist $2-\eta_r\sim C\,r^{-3/2}$ mit einer Log-Korrektur?' → 'Teilweise beantwortet: $(2-\eta_r)\,r\ln r$ ist fast konstant ≈ 2.20 (kn.verm.eta_r_log); dann geht der lokale Exponent wie $1+1/\ln r$ gegen 1. Prüfstein η₉.'
2. **`kn.beob.tpid10_gesetz`**: Kante ergänzen `{typ: verwandt, ziel: kn.verm.eta_r_log}`; in `bewertung.richtigkeit` anhängen: '; mit sechs Gliedern passt C/(r ln r) besser (kn.verm.eta_r_log)'.
3. **`kn.satz.kritische_rate`**: `pruefung.f5_neue_frage` ergänzen: 'Nächstes Glied: ι_r ln k (kn.satz.kritische_rate_log)'; in `aussage` nach β₈ einfügen: '$\beta_4=0.29348$, $\beta_5=0.37699$, $\beta_6=0.44293$ (kn.beob.beta_raenge)'.

## Strang und Ideenkarte

- **Strang `kn.etas`** (bestehend): Mitglieder + kn.satz.kritische_rate_log, kn.beob.iota_raenge, kn.beob.beta_raenge, kn.verm.eta_r_log, kn.gescheitert.e_potenzgesetz. `verlauf` anhängen: „Jede Berührung hat Krümmung β_r und iteratives Residuum ι_r (ab r = 4 negativ, um −0.6); 2 − η_r passt besser zu C/(r ln r) als zu einem Potenzgesetz — η₉ würde entscheiden; für e_r passt keine einfache Form.“
- **Ideenkarte**: bestehende Unteridee `stufen.etas` (Strang kn.etas) — keine neue Idee nötig; `kurz` optional: „Die Grenzbasen der höheren Stufen, ihr Grenzwert (TPID 10) und die Form jeder Berührung.“

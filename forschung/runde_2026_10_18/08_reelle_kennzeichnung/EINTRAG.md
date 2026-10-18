# Eintrag in den Graphen — Bericht 08 „Lokale, reelle Kennzeichnung“

Dateien: `daten/aussagen/r18_08_br.yaml`, `daten/theorie/r18_08_br.yaml` (9 Aussagen, je ein Theorie-Block).

## Angelegte ids
- br.satz.taylorradius_halbebene | ρ_T(x) ≤ x+2; Gleichheit überall ⇔ holomorph auf Re z > −2 ⇔ hohe Ableitungen wie log(x+2); Kneser maximal (A, B, E; bewiesen)
- br.satz.ungerade_ableitungen_halbebene | ungerade Ableitungen ≥ 0 ab Ordnung n₀ ⇒ Halbebenen-Holomorphie (C; bewiesen)
- br.satz.lokal_genau_dann_halbebene | lokale Kennzeichnung existiert ⇔ Q* (D; bewiesen)
- br.frage.halbebene_eindeutigkeit | Q*: erzwingt Holomorphie auf Re z > −2 schon Kneser? (offen; mit Teilfrage sup ρ_T < ∞)
- br.satz.slog_nie_total_monoton | slog′ nie total monoton, nie log-konkav/log-konvex (N1, N2; bewiesen)
- br.satz.endliche_vorzeichen_offen | Ränder universell; endlich viele Vorzeichen trennen nie (N3, N4; Skizze, Lücke Einfachheit der Nullstellen)
- br.beob.verbiegungen_singularitaeten | 13 Verbiegungen: Verzweigungspunkte in H, Radius 0.2–1.7, Vorzeichenbruch ab Ordnung 3–15 (E1–E3, X; numerisch)
- br.beob.kneser_ungerade_99 | Kneser: ungerade Koeffizienten bis Ordnung 99 positiv auf [−1.95, 1.5] (K; numerisch)
- br.beob.kneser_logkonvex | Kneser T′ log-konvex, min 0.9196; nicht kennzeichnend (L; numerisch)

Hinweis: gr.verm.ungerade_ableitungen (fremde Datei, nicht geändert) könnte im Text auf Ordnung 99 aktualisiert werden; die neue Messung steht als br.beob.kneser_ungerade_99 mit Kante `stuetzt`.

## Vorschlag Stränge
- Bestehender Strang `br.reell`: alle neun ids, Reihenfolge taylorradius → ungerade → lokal_genau_dann → frage → Negative (slog, endliche Vorzeichen) → Beobachtungen.
  verlauf: „18.10.: Lokale reelle Kennzeichnung auf eine Frage zurückgeführt — der Taylor-Radius x+2 an jeder Stelle kennzeichnet genau die Halbebenen-Holomorphie; ob diese Kneser erzwingt (Q*), ist offen.“
- Zusätzlich in `kn.eindeutigkeit`: br.frage.halbebene_eindeutigkeit (Paulsen–Cowgill ohne Grenzwert).

## Vorschlag Ideenkarte
- `fortsetzung.reell` (Reelle Kennzeichnung?): br.satz.taylorradius_halbebene und br.satz.lokal_genau_dann_halbebene als Kern; br.frage.halbebene_eindeutigkeit als offene Frage.
- `fortsetzung.eindeutigkeit`: br.frage.halbebene_eindeutigkeit.

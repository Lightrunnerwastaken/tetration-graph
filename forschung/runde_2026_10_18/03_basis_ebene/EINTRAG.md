# Eintrag Bericht 03 (Basis-Ebene) in den Graphen

Dateien: `daten/aussagen/r18_03_kl.yaml`, `r18_03_fp.yaml`, `r18_03_kn.yaml`; Theorie: `daten/theorie/r18_03_kl.yaml`, `r18_03_fp.yaml`, `r18_03_kn.yaml`.

## Angelegte ids
- kl.beob.kneser_fortsetzung_regulaer | Kneser über η = Fortsetzung der regulären Tetration (33 Stellen) (B3)
- kl.beob.eta_taylor | Taylor-Reihe von T_b(1/2) um η, Radius η − 1 (B4)
- kl.beob.regulaer_flach_singulaer | flache Singularität der regulären Familie bei λ = 1 (B5, B6, B9)
- kl.verm.eta_glatt_nicht_analytisch | C^∞, nicht analytisch; Defekt auf der regulären Seite (B7)
- kl.satz.c1_noerlund | erste Basisableitung bei η als Nörlund-Hauptlösung, Skizze (B8)
- kl.frage.eta_flacher_rest | ist Kneser exakt analytisch, Form von E⁻ (Offene Fragen 1–2)
- fp.satz.horn_eta | Kneser-Gebiet bei η ein Horn der Breite ∝ (b − η)^{3/2} (B1)
- kn.beob.motor_komplex_regulaer | Motor rechnet im Shell–Thron-Inneren regulär (B2)
- kn.beob.grosse_basis | T′(0), T″(0), Wendepunkt für b bis 10⁴ (B10)

## Vorschlag Stränge
- Bestehender Strang `kl.eta` („Der Grenzfall η“): aufnehmen kl.beob.kneser_fortsetzung_regulaer, kl.beob.eta_taylor, kl.beob.regulaer_flach_singulaer,
  kl.verm.eta_glatt_nicht_analytisch, kl.satz.c1_noerlund, kl.frage.eta_flacher_rest.
  Neue `verlauf`-Zeile: „Bei η passen regulär und Kneser zusammen: Kneser setzt die reguläre Reihe auf 33 Stellen fort; der Defekt ist eine flache Singularität auf der regulären Seite (glatt, nicht analytisch).“
- Bestehender Strang `fp.parameterraum`: aufnehmen fp.satz.horn_eta (verlauf unverändert).
- kn.beob.motor_komplex_regulaer: zu `kn.andere` oder dem Strang mit kn.satz.paulsen_komplex; kn.beob.grosse_basis: zu `gr.kurve` (Wendepunkte) — oder ohne Strang lassen.

## Vorschlag Ideenkarte
- Unteridee `fortsetzung.eta` („Der Grenzfall η“) trägt die kl-Einträge schon über `kl.eta`; `kurz` ergänzen: „… und warum die Tetration dort glatt, aber nicht analytisch in der Basis ist.“
- fp.satz.horn_eta über `turm.regionen` (Strang fp.parameterraum).

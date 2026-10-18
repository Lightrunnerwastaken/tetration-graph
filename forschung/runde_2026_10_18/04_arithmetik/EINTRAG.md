# Eintrag Bericht 04 (Arithmetik) in den Graphen

Dateien: `daten/aussagen/r18_04_ar.yaml`; Theorie: `daten/theorie/r18_04_ar.yaml`.

## Angelegte ids
- ar.satz.padisch_bewertung | exakte Bewertung V·(H − t) der p-adischen Turmziffern (B1, B2)
- ar.kor.padisch_hoehere_raenge | alle Ränge ≥ 4 mit gleichem p-adischem Grenzwert; D(3↑↑↑3) = 7 625 597 484 986 (B3, B4)
- ar.satz.padisch_lambert | p-adischer Turmgrenzwert als Lambert-Reihe in log_p a (B5)
- ar.satz.padisch_isometrie | a ↦ α_p(a) klassenweise analytische Isometrie (B6)
- ar.satz.padisch_koenigs | Multiplikator log_p α mit Bewertung V, Koenigs-Konstante κ_p(a), Ziffern-Gesetz (B7–B9)
- ar.satz.pentation_block_konstant | frischer Ziffernblock der Pentation friert ein; profiniter Turmgrenzwert α̂ (B10, B12)
- ar.satz.keine_padische_interpolation | keine p-adisch stetige Höhen-Fortsetzung ab Rang 4 (B11)
- ar.satz.ef_periode | Periode der exponentiellen Fakultät mod m ist Λ̃(m); Vorperiode (B13–B15)
- ar.konstr.ef_profinit | profinite exponentielle Fakultät, Ê(0) = 0, Ê(±1) = 1 (B16)

Nicht eingetragen: B17 (Stabilisierungshöhe ≈ Carmichael-Kettenlänge; numerisch, wohl bekannt, nicht verfolgt).
Quelle fehlt im Bestand: OEIS A127699 (im Text genannt, nicht als ref) — bei Bedarf in `daten/quellen` anlegen.

## Vorschlag Stränge
- Bestehender Strang `ar.kongruenzen` („Türme modulo k und p-adisch“): alle neun ids aufnehmen.
  Neue `verlauf`-Zeile: „Reste von Türmen werden konstant, p-adisch konvergiert jeder Turm — gegen eine Lambert-Reihe in log_p a, mit genau V·(H − t) fertigen Ziffern und einer Koenigs-Konstante für die nächsten; ab Rang 4 keine p-adische Höhen-Fortsetzung; die exponentielle Fakultät hat die Periode Λ̃(m).“
- Alternativ ar.satz.ef_periode und ar.konstr.ef_profinit in `ar.klassisch` neben ar.frage.ef_fortsetzung (dort stehen die EF-Einträge) — inhaltlich passt aber `ar.kongruenzen` besser.

## Vorschlag Ideenkarte
- Unteridee `zahlen.kongruenzen` (Straenge ar.kongruenzen, br.diskret); `kurz` ergänzen: „Türme modulo m, Endziffern, p-adische Grenzwerte als Lambert-Reihe und das Gesetz der nächsten Ziffern.“

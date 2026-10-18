# Eintrag Runde 18.10.2026 — 05 Kombinatorik

Dateien: `daten/aussagen/r18_05_kn.yaml`, `r18_05_it.yaml`, `r18_05_br.yaml`; Theorie `daten/theorie/r18_05_{kn,it,br}.yaml`.

## Angelegte ids
- kn.satz.stoerungsreihe_fluss | Störungsreihe = Fluss H′ = λᵤ(H) des Erzeugers von eˣ − 1 + u (S1, S10)
- kn.satz.stoerungsreihe_residuen | Lambert-Fixpunkte x± und Residuen 1/x± (S2)
- kn.satz.stoerungsreihe_bauform | H = γ + β tan Θ; Schichten polynomial in tan θ, ln cos θ, θ (S3, S4, S8)
- kn.satz.stoerungsreihe_logarithmen | höchste Logarithmen = 2 tan(θ + vL/3) (S5)
- kn.satz.stoerungsreihe_saekular | ω = v(1 − v⁴/270 + …), erstes nacktes θ in h₄ (S6)
- kn.satz.stoerungsreihe_dritte_schicht | h₂ geschlossen (S7)
- kn.verm.stoerungsreihe_ganzzahlig | normierte Schichtfolgen ganz (P3, Frage 3)
- it.satz.julia_variation_bernoulli | u-linearer Teil des Erzeugers über Bernoulli-Zahlen (S9)
- it.gescheitert.pentation_baeume | keine Baum-Kombinatorik für die Pentation (Q1)
- br.satz.produkt_iteration_ketten | Iterierte von x·eˣ = Ketten punktierter Partitionen, Möbius (−1)^{n−1}n^{n−2} (P1–P3)

Die Taylor-Befunde T1–T4 stehen gemeinsam mit Bericht 06 in `r18_06_kn.yaml` (kn.satz.minus2_hebbar, kn.satz.taylor_darboux, kn.beob.taylor_drei_stufen); dort sind beide Berichte zitiert.

## Vorschlag Stränge
- Neuer Strang `kn.stoerungsreihe` „Die Störungsreihe der Tetration“: mitglieder kn.konstr.stoerungsreihe, kn.satz.stoerungsreihe_tangens, kn.satz.stoerungsreihe_baeume, kn.satz.stoerungsreihe_fluss, kn.satz.stoerungsreihe_residuen, kn.satz.stoerungsreihe_bauform, kn.satz.stoerungsreihe_logarithmen, kn.satz.stoerungsreihe_saekular, kn.satz.stoerungsreihe_dritte_schicht, kn.verm.stoerungsreihe_ganzzahlig, it.satz.julia_variation_bernoulli, kn.frage.stoerungsreihe.
  verlauf: „Von der Rekursion über Bäume und Tangensschichten zur geschlossenen Bauform H = γ + β tan Θ (18.10.2026).“
- br.satz.produkt_iteration_ketten und it.gescheitert.pentation_baeume in den bestehenden Strang, der br.beob.produkt_iteration_dreieck / it.satz.zeta_polynom_partitionen enthält (Kombinatorik der Iteration).

## Vorschlag Ideenkarte
- Unter der Idee zur Störungsreihe bzw. zur Kombinatorik der Iteration (`daten/ideen.yaml`) eine Unteridee „Bauform der Störungsreihe“, kurz: „Alle Schichten aus einem verbogenen Tangens; offen bleibt der reguläre Rest R_u.“

## Fehlende Quellen (nicht angelegt)
Chapoton–Vallette, „Pointed and multi-pointed partitions of type A and B“ (arXiv math/0410051); Hoffman, Ableitungspolynome von tan/sec; Flajolet–Sedgewick, Analytic Combinatorics.

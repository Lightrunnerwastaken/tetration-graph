# Bericht Strang L2 (kn, bw) — Stand 12.10.2026, 15:45

## Erledigt
- Karten: KARTE_kn.md (Kneser, Verfahren, halbe exp, Hyperoperationen, TPID 10 mit allen Zugängen, ganze TPID-Liste mit Zuordnung, Lücken),
  KARTE_bw.md (Basiswechsel-Arbeiten anderer, Zugänge je Frage, Lücken).
- TPID 10: neuer Strang `kn.etas` mit 9 Aussagen (daten/aussagen/lit_kn_1.yaml): Interpolation (Robbins/mike3), Matrix-slog (nuninho), Levenstein 2017
  (η₅, η₆, η₇), Abstandsgesetz mit Prüfwert η₇, TPID 11, Nixons Beschränktheits-Heuristik, zwei Nebenfragen. kn.verm.tpid10 neu gegliedert.
- Experiment exp.tpid10_abstaende (e70, 8/8 Prüfungen, < 1 s, ohne Kneser-Motor).
- Neue Literaturknoten: kn.konstr.paulsen_stoerung (Paulsen 2016, Teile gelesen), bw.satz.paulsen_cowgill_hardy (Paulsen–Cowgill 2021, zbMATH).
- Zuschreibungen korrigiert (Befund auf dem BRETT): kn.beob.stufe5, stufe6, unendliche_hexation, hexation_struktur → bekannt; bw.lem.randgesetz → bekannt;
  Recherche dokumentiert bei sechs Hexations-Knoten und bw.satz.imaginaere_richtung (→ vermutlich_neu); Quellen ergänzt bei kn.konstr.pentation
  (Kouznetsov 2014), kn.beob.tetra_euler, kn.frage.hexation_kneser, bw.satz.phasengesetz, bw.verm.nirgends_analytisch, bw.frage.halbiter_vergleich.
- `--nur-pruefen`: 0 Fehler, 0 Warnungen (auch die alte Warnung bei bw.kor.modenzahl durch `nachgeprueft` erledigt).

## Entscheidungen
- Forumsbeiträge, die über den Seitenabruf gelesen wurden: geprueft original; nur Suchtreffer: sekundaer.
- Paulsen 2016: Seiten 81–83 und 102–105 gelesen (per-Zitat `geprueft: original`); Paulsen–Cowgill 2021 und Paulsen 2026 nicht zugänglich.
- Neue Etas-Knoten in eigenem Strang `kn.etas` (Lehrbuch-Platz bei L3 angefragt); kn.beob.tetra_euler, stufe5, stufe6 bleiben in kn.hoehere_stufen.

## Offene Lücken (Priorität A)
- Paulsen 2026 und Paulsen–Cowgill 2021 im Volltext (Leitung fragen).
- η₇, e₇ genau rechnen (Kneser-Motor, Heptation über die Kneser-artige Hexation) — prüft Abstandsgesetz und TPID 11; e₇: Hochrechnung 3.54 gegen Levensteins „nahe 3.7“.
- Matrix-slog bei η₄ mit wachsender Dimension (kn.beob.eta4_matrix f5).

## Block 2 (12.10.2026, 16:00)
- Neu (10): kn.satz.tk_kriterien, kn.verm.tk_c2, kn.konstr.paulsen_cowgill_iteration, kn.beob.kouznetsov_numerik, kn.satz.walker_konjugation,
  kn.verm.szekeres_wachstum, kn.verm.tpid3, kn.verm.tpid3_bedingung, kn.frage.tpid14 (lit_kn_2.yaml), bw.frage.nixon_eta_basiswechsel (lit_bw_1.yaml).
- Verbessert (14): kn.satz.trappmann_kouznetsov (Satz 1/2 nach dem Original), kn.satz.paulsen_cowgill, kn.satz.paulsen_komplex (Holomorphie in b nicht belegt),
  kn.satz.kneser_eigenschaften, kn.konstr.kneser, kn.konstr.walker, kn.konstr.carleman (Walker vor Robbins), kn.konstr.cauchy, kn.satz.singulaer_minus2 und
  kn.satz.linear_glatt (folklore → bekannt), kn.konstr.linear, kn.frage.beta_methode, kn.verm.carleman_eigenwerte, kn.beob.referenzwerte.
- Experiment exp.tpid3 (e71, 5/5, 63 s). PDFs (Paulsen 2016, Trappmann–Kouznetsov 2011) legt WebFetch selbsttätig im Sitzungsordner ab, nicht im Projekt.

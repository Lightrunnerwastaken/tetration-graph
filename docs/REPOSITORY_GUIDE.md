# Repository guide

The knowledge graph is written in German and kept in sync in English. File, folder and field names are therefore German.
This guide explains what everything means, so you can read, check and extend the graph without knowing German.

## Layout

| Path | Meaning | Contents |
|---|---|---|
| `baue.py` | *build* | checks all data and builds `ausgabe/tetration-graph.html`; `--nur-pruefen` = check only, `--experimente` = recompute experiments first |
| `daten/` | *data* | everything the graph consists of (YAML) |
| `daten/aussagen/*.yaml` | *statements* | one file per area (`00_grundlagen` … `10_…`) plus supplements (`lit_*`, `kn_*`) |
| `daten/theorie/*.yaml` | *theory* | the textbook explanation of each statement (parts `worum` = what it is about, `genau` = precisely, `warum` = why it holds, `beispiel` = example, `einordnung` = context) |
| `daten/quellen.yaml`, `daten/quellen/` | *sources* | literature, forum threads, own papers |
| `daten/experimente.yaml`, `daten/experimente/` | *experiments* | registry: script, result file, figure, supported statements |
| `daten/objekte.yaml` | *objects* | vocabulary of mathematical objects (fixed point L, Kneser tetration, …) |
| `daten/straenge/` | *strands* | research strands per area (reading paths) |
| `daten/ideen.yaml` | *ideas* | the idea map: 8 main ideas, their sub-ideas, and which strands belong where |
| `daten/lehrbuch.yaml` | *textbook* | chapters and their order |
| `daten/gebiete.yaml` | *areas* | the 9 areas (see below) |
| `daten/uebersetzung/`, `daten/uebersetzung_struktur/` | *translation* | English versions of statements and of all structure texts, with a hash of the German source |
| `wissen/` | *knowledge* | Python package: `laden` = loading, `graph` = checking and status propagation, `schema` = vocabulary, `landkarte` = idea map, `lehrbuch` = textbook, `uebersetzung` = translation tools, `ausgabe` = output |
| `experimente/e*.py` | *experiments* | reproducible computations with self-checks; results in `ergebnisse/`, figures in `bilder/` |
| `ansicht/` | *view* | HTML template, style and JavaScript modules of the viewer |
| `ausgabe/` | *output* | built files: `tetration-graph.html` (the viewer), `pruefbericht.md` (audit report), `landkarte.md` (idea map as text), `graph.json` |
| `forschung/` | *research* | working folders of research rounds; only established results are taken over into the graph |
| `docs/` | | reports, design documents, translation glossary (`docs/uebersetzung/`), development log (`HISTORY.md`) |
| `tests/` | | tests of the checker |

### Areas (`gebiet`)

| Prefix | German | English |
|---|---|---|
| `gr` | Grundlagen | foundations |
| `fg` | Funktionalgleichungen | functional equations |
| `it` | Lokale Iterationstheorie | local iteration theory |
| `fp` | Fixpunkte und unendlicher Turm | fixed points and the infinite tower |
| `kl` | Lösungsklassen nach Basis | classes of solutions by base |
| `kn` | Kneser-Tetration und slog | Kneser's construction, super-logarithm, higher hyperoperations |
| `bw` | Basiswechsel | base change |
| `ar` | Arithmetik und Transzendenz | arithmetic and transcendence |
| `br` | Brücken zu anderen Gebieten | bridges to other fields |

A statement id is `<area>.<kind>.<name>`, e.g. `kn.satz.kneser_eigenschaften`.

## Fields of a statement

| Field | Meaning |
|---|---|
| `id`, `titel` | identifier, title |
| `art` | kind: `definition`, `satz` theorem, `lemma`, `korollar` corollary, `konstruktion` construction, `annahme` assumption, `vermutung` conjecture, `beobachtung` observation, `gegenbeispiel` counterexample, `gescheitert` failed approach, `korrektur` correction, `frage` open question |
| `ebene` | layer: **A** every solution of the functional equation · **B** a class of solutions · **C** a concrete construction · **X** auxiliary result from another field |
| `status` | see the next table |
| `niveau` | level: `einstieg` entry (school mathematics + calculus), `vertiefung` university, `forschung` research |
| `objekte` | objects involved (keys of `daten/objekte.yaml`) |
| `bereich` | domain of validity |
| `voraussetzungen` | hypotheses, each with `benutzt: ja/nein/unklar` (used in the proof: yes/no/unclear) |
| `aussage` | the statement (Markdown + TeX) |
| `einfach` | the same in plain words |
| `formen` | equivalent formulas (`tex`, `notiz` = note) |
| `beweis` | proof: `art` (`vollstaendig` complete, `literatur`, `skizze` sketch, `bedingt` conditional), `text`, `luecken` = gaps |
| `stuetzt_sich_auf` | the proof relies on these statements |
| `bedingungen` | statements that are assumed (for conditional results) |
| `kanten` | further edges (see below) |
| `quellen` | sources: `ref` (key in `daten/quellen*`), `stelle` = location in the source |
| `experimente` | experiments supporting it |
| `pruefung` | five check questions: `f1_folgt_unmittelbar` does it follow immediately? · `f2_benutzte_voraussetzungen` which hypotheses are really used? · `f3_spezialfall` special case of something more general? · `f4_uebertragbar` transferable to other constructions? · `f5_neue_frage` what new question does it raise? |
| `bewertung` | assessment, kept separate: `richtigkeit` correctness, `neuheit` novelty, `nutzen` usefulness (+ `nutzen_stufe` hoch/mittel/gering = high/medium/low), `recherche` = documented literature search |
| `naechste_frage`, `gegenbeispielsuche` | next question; documented search for counterexamples (required for conjectures) |

### Status

| `status` | Meaning |
|---|---|
| `definition` | definition |
| `bewiesen_literatur` | proved in the literature |
| `bewiesen` | proved here / in own work |
| `entlastet` | assumption discharged (now proved) |
| `bedingt` | proved under conditions that are not yet proved |
| `skizze` | proof sketch with named gaps |
| `numerisch` | observed numerically |
| `vermutung` | conjecture |
| `annahme` | unproved assumption |
| `offen` | open |
| `gescheitert` | failed approach |
| `widerlegt` | refuted |
| `ueberholt` | superseded |

The checker computes an **effective status**: a `bewiesen` statement that relies on something weaker becomes `bedingt`, and the unproved ends of the chain are listed.

### Novelty (`neuheit`) and source check (`geprueft`)

`klassisch` classical · `bekannt` known from the literature · `folklore` probably known, no source found · `vermutlich_neu` probably new after a documented search ·
`neu_ungeprueft` called new, search not documented · `nur_schreibweise` new notation only · `unklar` unclear.
Sources are marked `original` (checked against the original), `sekundaer` (via a summary or review) or `nein` (not checked).

### Edge types (`kanten.typ`)

`stuetzt_sich_auf` relies on (proof) · `folgt_unmittelbar_aus` follows immediately from · `spezialfall_von` special case of · `verallgemeinert` generalizes ·
`ersetzt` supersedes · `korrigiert` corrects · `widerlegt` refutes · `entlastet` discharges an assumption · `stuetzt` supports (evidence) ·
`uebertragbar_auf` transferable to · `analog_zu` analogous to · `motiviert` motivates · `verwandt` related.

## What the checker enforces

- references and vocabulary exist; no cycles in proof chains;
- status propagation (proved on a weaker basis → conditional; root conditions named);
- anything depending on a refuted or superseded statement is flagged “re-check” (`neu-pruefen`);
- conjectures need a documented counterexample search, sketches named gaps, “probably new” a documented literature search;
- an experiment result counts only if the hash of its script matches;
- every statement belongs to a strand, a textbook chapter and an idea; translations are up to date and keep formulas, numbers and references unchanged.

Findings are written to `ausgabe/pruefbericht.md` as errors (`FEHLER`), warnings (`WARNUNG`) and hints (`HINWEIS`).

## Adding a statement

1. Add it to a file in `daten/aussagen/` (template below), with a theory block in `daten/theorie/`.
2. Put its id into a strand (`daten/straenge/<area>.yaml`); new strands also need a place in `daten/lehrbuch.yaml` and `daten/ideen.yaml`.
3. Run `python baue.py --nur-pruefen` until there are no errors or warnings.
4. Add the English version: `python -m wissen.uebersetzung export --gebiet <area> --felder kern --nur-fehlend -o batch.json`, translate, then
   `python -m wissen.uebersetzung einlesen <file> --modell <name>` (details and glossary: `docs/uebersetzung/`).

```yaml
- id: kn.satz.example
  titel: Short title
  art: satz
  ebene: B
  status: bewiesen
  objekte: [kneser_T]
  bereich: '$b>\eta$'
  voraussetzungen:
    - {text: 'T strictly increasing', benutzt: ja}
  aussage: |
    Statement with $\TeX$.
  einfach: >-
    The same in plain words.
  beweis: {art: vollstaendig, text: '…', luecken: []}
  stuetzt_sich_auf: [kn.konstr.kneser]
  quellen: [{ref: kneser1950, stelle: 'Theorem 1'}]
  pruefung: {f1_folgt_unmittelbar: '…', f2_benutzte_voraussetzungen: '…', f3_spezialfall: '…', f4_uebertragbar: '…', f5_neue_frage: '…'}
  bewertung: {richtigkeit: '…', neuheit: bekannt, nutzen: '…', nutzen_stufe: mittel}
```

(The German texts are the source; please write new statements in German if you can, or in English — they will be translated.)

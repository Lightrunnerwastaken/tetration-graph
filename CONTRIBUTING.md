# Contributing

Thank you for helping to make this map of tetration more complete and more correct. Three kinds of contributions are especially useful.

## 1. Report an error

If a statement, a proof, a number or a citation is wrong, please [open an issue](https://github.com/Lightrunnerwastaken/tetration-graph/issues) with

- the **statement id** (shown in the detail panel of the viewer, e.g. `kn.satz.kneser_eigenschaften`),
- what is wrong, and
- if possible, a source or a short calculation.

Corrections are recorded, not hidden: the old version stays visible as superseded or refuted, and everything that depended on it is flagged for re-checking.

## 2. Point to missing literature

If a result in the graph is already known — or a known result is missing — please open an issue with the reference (author, title, journal or arXiv/forum link)
and, ideally, the exact place (theorem number, page, post). Attribution matters to us: several results here were re-attributed after such hints.

## 3. Propose a statement or a correction directly

Pull requests are welcome. The data are YAML files; [docs/REPOSITORY_GUIDE.md](docs/REPOSITORY_GUIDE.md) explains every field (the names are German) and contains a template.
Before opening a pull request:

```bash
pip install pyyaml pytest
python baue.py --nur-pruefen     # must report 0 errors and 0 warnings
python -m pytest -q tests
```

Please state honestly what is proved, what is only numerical, and what is conjectured — the status field is the heart of the graph.
You may write in English or German; the project keeps a German source text and an English translation, and we will take care of the other language.

## Ground rules

- Every claim needs a source, a proof, or a reproducible experiment.
- Only cite what you have actually read; mark second-hand references as such.
- Be kind. Mathematics is hard, and so is checking it.

---

**Deutsch, kurz:** Fehler und fehlende Literatur gern als Issue melden (mit der Kennung der Aussage). Pull Requests sind willkommen; vorher `python baue.py --nur-pruefen` (0 Fehler, 0 Warnungen) und die Tests laufen lassen. Aufbau und Feldnamen: [docs/REPOSITORY_GUIDE.md](docs/REPOSITORY_GUIDE.md).

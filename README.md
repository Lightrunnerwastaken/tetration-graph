# Tetration Knowledge Graph

**A verifiable map of what is known about tetration** — the operation $a\uparrow\uparrow x$, “$a^{a^{\cdots}}$ with $x$ storeys”, extended to real and complex heights —
and of its connections to functional equations, iteration theory, complex dynamics, arithmetic and logic.

**[Open the interactive graph](https://lightrunnerwastaken.github.io/tetration-graph/)** · runs in the browser, no installation · Deutsch: [README.de.md](README.de.md)

[![Idea map of the knowledge graph](docs/bilder/ideen_karte.png)](https://lightrunnerwastaken.github.io/tetration-graph/)

Every entry is a single mathematical **statement** with its hypotheses, proof or proof gaps, sources, and an explicit status
(proved · proved in the literature · conditional · numerical · conjecture · open · refuted …). A checker propagates these statuses along proof chains:
a theorem that rests on an unproved assumption is automatically marked *conditional*, and anything depending on a refuted result is flagged for re-checking.
Numerical claims are backed by reproducible experiments whose results are tied to the hash of the script that produced them.

**Current size:** about 390 statements, 680 edges, 55 experiments, 300 sources, in 9 areas and 8 main ideas.
Open the graph for the authoritative numbers, or see the audit report [`ausgabe/pruefbericht.md`](ausgabe/pruefbericht.md).

## What you can do with it

- **Ideas** (start view): the structure of the field as a network of eight main ideas — the functional equation, dynamics of $b^z$, the infinite tower,
  fractional continuation, base change, higher hyperoperations, arithmetic, bridges to logic and analysis. Their position is computed from how strongly
  their statements depend on each other; rings show how much is proved and where it comes from (literature, forum, own work). Click to zoom in to sub-ideas and statements.
- **Map:** a zoomable map by area and research strand.
- **Theory:** a textbook in 12 chapters, from first examples to current research, with a worked explanation for every statement.
- **Research:** open questions, gaps, corrections and experiments.
- For conjectures and conditional theorems, the detail panel shows **what the result is made of** — which pieces are proved and where exactly something is missing.

## Some results in the graph

- **Literature map:** results of Kneser, Écalle, Paulsen–Cowgill, Trappmann–Kouznetsov, Levenstein and many others are entered with attribution;
  several statements first believed new were traced to earlier work and marked as known (e.g. a series of critical bases already computed by Levenstein in 2017).
- **Critical bases of higher hyperoperations:** $\eta_7=1.8397586607\ldots$ and $\eta_8=1.8681465773\ldots$, extending the sequence
  $e^{1/e},1.6353,1.7374,1.7991,\ldots$ that is conjectured to tend to 2 (TPID 10) — together with the proof that the sequence depends on the chosen continuation.
- **Polynomial square roots of exp:** polynomials $p_N$ with $p_N(p_N(x))=e^x+O(x^{N+1})$ approach Kneser's half-exponential function;
  proved conditionally, with the remaining core problem isolated.
- **Half rank:** first values of a hyperoperation “3.5” between $2^x$ and tetration at base 2 (about 4 digits), and the reason why
  perturbation theory cannot give more: the remaining freedom is exactly the periodic θ-freedom of every fractional iteration.
- **Mixed-base tetration** (base change $T_b\circ\mathrm{slog}_d$): phase law, circle map, Fourier modes and certified reconstruction.

## How it is organized

The data is plain YAML; the build script checks everything and produces a single offline HTML file.
**The project is written in German and translated into English**, so file, folder and field names are German
(`daten/aussagen` = statements, `beweis` = proof, `quellen` = sources, …).
**[docs/REPOSITORY_GUIDE.md](docs/REPOSITORY_GUIDE.md)** explains the layout, every German name, the status vocabulary and how to add a statement.

```bash
pip install pyyaml pytest          # enough to check and build
python baue.py                     # check + build → ausgabe/tetration-graph.html
python baue.py --nur-pruefen       # audit report only → ausgabe/pruefbericht.md
python -m pytest -q tests          # tests of the checker
```

The experiments additionally need `mpmath`, `numpy`, `scipy`, `matplotlib`; some (the Kneser engine, hexation) also PARI/GP.

## Contributing and citing

Corrections, missing literature and new results are welcome — see **[CONTRIBUTING.md](CONTRIBUTING.md)**.
To cite the graph, use the “Cite this repository” button on GitHub (from [CITATION.cff](CITATION.cff)).
The development log of the project is in [docs/HISTORY.md](docs/HISTORY.md).

## License

Code under the MIT license ([LICENSE](LICENSE)); content (statements, texts, figures) under CC BY 4.0 ([LICENSE-INHALT.md](LICENSE-INHALT.md)).

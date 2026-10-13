# Literaturkarte kn — Kneser-Tetration, Superlogarithmus, Rechenverfahren, höhere Hyperoperationen

Strang L2, Stand 12.10.2026. Suchwege: Bestand (05_kneser_slog.yaml, quellen.yaml), Tetration-Forum (TPID-Liste tid 162 alle vier Seiten,
Bibliografie-Thread „Major references“ tid 1681, Themen zu η₄/e₄, Pentation, Hexation), Hyperops-Wiki, zbMATH Open (API, „tetration“),
arXiv (API, „tetration“, „pentation/hyperoperation“), Crossref, OEIS (Ziffernsuche), MathOverflow/Math.SE (Websuche), Zitationsketten
(Forumsbeiträge, die Arbeiten nennen).

**Prüfstand** (Spalte P): **O** = Original (oder die zitierte Stelle) gelesen · **S** = über Abstract, Rezension, Metadaten oder zitierenden
Beitrag geprüft · **N** = nur aus Bibliografie-Liste oder Gedächtnis, vor Verwendung nachschlagen.
**Bei uns:** Knoten-ID bzw. Quellenschlüssel; „—“ = fehlt.

## 1. Knesers Konstruktion, Eindeutigkeit, komplexe Basen

| Autor(en), Jahr | Werk, Ort | Hauptergebnis | Methode | Beweisstand | P | Bei uns |
|---|---|---|---|---|---|---|
| H. Kneser 1950 | Reelle analytische Lösungen der Gleichung φ(φ(x)) = eˣ …, J. reine angew. Math. 187, 56–67 | reell-analytische halbe Iterierte von eˣ, damit reell-analytische Tetration zur Basis e | Schröder-Funktion am komplexen Fixpunkt L, Riemann-Abbildung korrigiert zu reell | bewiesen | S | kn.konstr.kneser, kn.satz.halbe_exp, kneser1950 |
| D. Kouznetsov 2009 | Solution of F(z+1)=exp(F(z)) in complex z-plane, Math. Comp. 78, 1647–1670 | Superexponential zur Basis e in der ganzen Ebene ohne (−∞,−2], Werte, Bilder | Cauchy-Integral im Streifen, Iteration | numerisch (Existenz nicht bewiesen) | S (zbMATH) | kn.konstr.cauchy, kn.beob.kouznetsov_numerik |
| D. Kouznetsov 2010 | Tetrational as special function, Vladikavkaz Mat. Zh. 12(2), 31–45 | Tetration als spezielle Funktion: Eigenschaften, Asymptotik, Approximationen | wie 2009, Darstellung als Spezialfunktion | numerisch/Übersicht | S (zbMATH) | Quelle bei kn.konstr.cauchy |
| H. Trappmann, D. Kouznetsov 2011 | Uniqueness of holomorphic Abel functions at a complex fixed point pair, Aequationes Math. 81, 65–76 | Eindeutigkeit der Abel-Funktion mit Grenzverhalten an L, L̄ | Funktionentheorie (Streifen, periodische Störung) | bewiesen | O (arXiv:1006.3981, S. 3–9, 14) | kn.satz.trappmann_kouznetsov (Satz 1/2 genau), kn.satz.tk_kriterien, kn.verm.tk_c2 (Vermutung S. 4) |
| D. Kouznetsov, H. Trappmann 2010 | Portrait of the four regular super-exponentials to base √2, Math. Comp. 79, 1727–1756 | vier reguläre Superexponentiale unter η | reguläre Iteration | bewiesen/numerisch | S | kouznetsov_trappmann2010 (Gebiet kl, L1) |
| H. Trappmann, D. Kouznetsov 2012 | Computation of the two regular super-exponentials to base e^{1/e}, Math. Comp. 81, 2207–2227 | Basis η, parabolisch | Abel-Funktion am parabolischen Fixpunkt | bewiesen/numerisch | S | kouznetsov_trappmann2012 (L1) |
| W. Paulsen, S. Cowgill 2017 | Solving F(z+1)=b^{F(z)} in the complex plane, Adv. Comput. Math. 43, 1261–1282 | Eindeutigkeit: holomorph in der rechten/oberen Halbebene mit Grenzwert L | Funktionentheorie + Iterationsverfahren | bewiesen | S (Abstract; Nixon 2021 tid 162 #21) | kn.satz.paulsen_cowgill, kn.konstr.paulsen_cowgill_iteration (180 Knoten, > 50 Stellen) |
| S. Cowgill 2017 | Exploring tetration in the complex plane, Masterarbeit, Arkansas State Univ. | Vorstufe zu 2017/2019 | — | — | N | — |
| W. Paulsen 2019 | Tetration for complex bases, Adv. Comput. Math. 45, 243–267 | Eindeutigkeit und Berechnung für komplexe Basen, auch innerhalb Shell–Thron | Iteration | bewiesen/numerisch | S (Abstract) | kn.satz.paulsen_komplex — Holomorphie in b dort **nicht belegt** |
| W. Paulsen 2026 | Analyzing the function z↑↑a for fractional a, Adv. Comput. Math. 52(2), Nr. 25 (online 01.04.2026) | Tetration fester gebrochener Höhe als Funktion der **Basis** | — | — | S (nur Metadaten, Crossref/zbMATH; Inhalt nicht zugänglich) | — (Lücke A: wichtig für bw.ann.basis_diffbar) |
| S. Levenstein 2010–2011 | kneser.gp (2010), fatou.gp (2011), Tetration Forum tid 1017 | Kneser auf 30–100 Stellen, auch komplexe Basen | θ-Iteration auf einem Kreis, Matrix-Lösung | numerisch | O | levenstein_fatou, gp_tetration |
| V. Reshetnikov 2017 | MathOverflow 259278 | explizite Reihe für analytische Tetration (Basen ≤ η) | Reihen | numerisch | O | reshetnikov2017 |
| J. D. Nixon 2021 (Forum 23.08.2021, tid 162) | Eindeutigkeitsfrage im Stil Paulsen–Cowgill (ohne TPID-Nummer) | offen | — | Frage | S | — |

## 2. Andere Konstruktionen und Rechenverfahren

| Autor(en), Jahr | Werk, Ort | Hauptergebnis | Methode | Beweisstand | P | Bei uns |
|---|---|---|---|---|---|---|
| C. W. Clenshaw, D. W. Lozier, F. W. J. Olver, P. R. Turner 1986 | Generalized exponential and logarithmic functions, Comput. Math. Appl. 12B | verallgemeinerte exp/log für Level-Index-Arithmetik (stückweise, „linear“) | Definition über Iteration + Intervall | bewiesen | S (zbMATH) | Quelle bei kn.konstr.linear |
| P. L. Walker 1991a | Infinitely differentiable generalized logarithmic and exponential functions, Math. Comp. 57, 723–733 | C^∞-Tetration/slog zur Basis e | Grenzwert dexp^{−n}∘expⁿ (dexp = eˣ−1) | C^∞ bewiesen, Analytizität offen | S (zbMATH) | kn.konstr.walker, kn.satz.walker_konjugation |
| P. L. Walker 1991b | On the solutions of an Abelian functional equation, J. Math. Anal. Appl. 155, 93–110 | zwei Lösungen von f(e^{z/e}) = f(z)+1 (Basis η) | Iteration von 1−e^{−z}, log(1+z) | bewiesen | S (zbMATH) | walker1991_jmaa; nach Trappmann–Kouznetsov 2011 auch Matrix-Koeffizienten des 4-Logarithmus (vor Robbins) → kn.konstr.carleman |
| G. Szekeres 1962 (Jahr nach zbMATH, Hinweis L1) | Fractional iteration of exponentially growing functions, J. Austral. Math. Soc. 2, 301–320 | Abel-Funktion für exp-artige Funktionen, „reguläres Wachstum“ | Asymptotik bei ∞ | bewiesen | S | szekeres1961 (L1: it) |
| G. Szekeres 1998 | Abel's equation and regular growth: variations on a theme by Abel, Experiment. Math. 7(2), 85–100 | Abel-Funktion von eˣ und Wachstumsskala, Numerik | Iteration von eˣ−1 / Asymptotik | teils numerisch | S (zbMATH) | kn.verm.szekeres_wachstum |
| A. Robbins 2005 | Solving for the analytic piecewise extension of tetration and the super-logarithm (selbst veröffentlicht) | „natürlicher“ slog | abgeschnittenes lineares System (Carleman/Abel-Matrix) | numerisch | N | kn.konstr.carleman, robbins2005 |
| H. Trappmann 2010 | The intuitive logarithm, arXiv:1008.1409 | intuitive Abel-Funktion von bx ist log_b (TPID 9 gelöst) | Abel-Matrix | bewiesen | S | kn.satz.intuitiver_log |
| G. Helms 2007– | Matrix-Methoden, Eigenzerlegung, „tetdocs“ (go.helms-net.de/math/tetdocs), Forum | Potenz der Carleman-Matrix, Basiswechsel unter η, TPID 18 (2022) | Eigenzerlegung abgeschnittener Matrizen | numerisch | S | tetrationforum, kn.verm.carleman_kneser |
| R. Aldrovandi, L. P. Freitas 1998 | Continuous iteration of dynamical maps, J. Math. Phys. 39 | Bell-Matrix-Iteration | Bell-/Carleman-Matrizen | formal | N | — |
| R. Aldrovandi 2014 | Tetration: an iterative approach, arXiv:1410.3896 | Tetration über Bell-Matrizen | Bell-Matrizen | formal/numerisch | S (arXiv-Liste) | — |
| J. Fox 2007 | Forum (Matrix-slog, „Accelerated Abel matrix“) | Robbins-slog mit 1200 Gliedern; Pentation e↑↑↑x | Matrix-slog | numerisch | O | forum_fox2007_pentation |
| Ansus, mike3 2009–2010 | Forum tid 343, 353, 437 | Kontinuumssumme/-produkt | Kontinuumssumme | numerisch | O | kn.konstr.kontinuumssumme |
| tommy1729 2010–2017 | Forum tid 424, 1144 | 2sinh-Methode | Grenzwert ln^[n]∘2sinh^[n+x] | numerisch | O | kn.konstr.zwei_sinh |
| J. D. Nixon 2021–2022 | arXiv 2104.01990, 2101.03021, 2108.13519, 2208.05328; Forum tid 1334 | β-Methode: asymptotische Lösungen über unendliche Kompositionen | unendliche Kompositionen | teils bewiesen, Verhältnis zu Kneser umstritten | S | kn.frage.beta_methode, nixon2021b |
| MphLee, Nixon 2021 | Forum tid 1290 | verallgemeinerter Kneser-Grenzwert-Trick | iterierte Grenzwerte | Skizze | O | kn.satz.grenzwert_trick |
| W. Paulsen 2016 | Finding the natural solution to f(f(x)) = exp(x), Korean J. Math. 24, 81–106 | Störungsmethode: Abel-Funktion einer gestörten Exponentialfunktion mit reellem Fixpunkt, Grenzübergang; Metrik zwischen Lösungen; neue Lösung „zwischen“ den bekannten; Stetigkeit bei a → η (Kor. 1); Engpass-Zählung π√2/√(η_a f''(x_a)) im Beweis von Lemma 2 | Störung + Grenzwert | bewiesen | O (S. 81–83, 102–105) | paulsen2016; neu: kn.konstr.paulsen_stoerung |
| T. Ueda 2021 | Extension of tetration to real and complex heights, arXiv:2105.00247 | weitere Fortsetzung reeller/komplexer Höhen | — | — | S (arXiv-Liste) | — |
| S. Nesargi, G. Roudenko 2025 | Analytical and numerical approaches for finding functional iterates and roots, arXiv:2509.24049 | slog und Superwurzel als Werkzeug für gebrochene Iterierte von eˣ, Vergleich mit Kneser | numerisch | numerisch | S (Suchtreffer) | — |
| R. Munafo ~2004 | Webseite hyper4 (mrob.com/pub/math/hyper4.html) | eigene Fortsetzung auf reelle Höhen (Potenzreihe in der Höhe) | Turm-Methode | heuristisch | S (Robbins' Übersicht) | — |
| I. Galidakis ~2004 | Webseiten zu hyper4 und unendlichen Exponentialen, Literaturliste IERefs | stückweise Fortsetzungen; Bibliografie | — | — | N | — |
| H. Trappmann, A. Robbins 2007–2008 | Tetration FAQ (Fassungen 09.08.2007, 12.01.2008, 10.07.2008), Forum tid 184 | Übersicht aller Fortsetzungen | Übersicht | — | S | — (Lücke C: als Sekundärquelle nutzen) |
| H. Trappmann, D. Kouznetsov 2009/2010 | „5+ methods for real analytic tetration“ (Vorabdruck, im Bibliografie-Thread tid 1681) | Vergleich von fünf reell-analytischen Verfahren | Vergleich | numerisch | N (nicht auffindbar) | — |
| K. Kefalas 2014 | On smooth solutions of non-linear dynamical systems f_{n+1}=u(f_n), Part I | glatte Lösungen, zitiert Hooshmand, Walker, Szekeres | — | — | N | — |
| M. H. Hooshmand 2006 | Ultra power and ultra exponential functions, Integral Transforms Spec. Funct. 17 | eindeutige „ultra exponential“ unter Konvexitätsbedingung (linear auf [−1,0]) | Charakterisierung | bewiesen | S | kn.satz.hooshmand |

## 3. Halbe Exponentialfunktion und Wachstum

| Autor(en), Jahr | Werk, Ort | Hauptergebnis | Methode | P | Bei uns |
|---|---|---|---|---|---|
| G. H. Hardy 1910 | Orders of Infinity | Wachstumsskalen | Hardy-Körper | S | hardy1910 |
| M. Boshernitzan 1986 | Hardy fields and existence of transexponential functions, Aequationes Math. 30 | transexponentielle Funktionen in Hardy-Körpern | Hardy-Körper | S | boshernitzan1986 |
| J. van der Hoeven 2006; G. A. Edgar 2011 | Transseries …, LNM 1888; Fractional iteration of series and transseries | keine halbe exp in Transreihen | Transreihen | S / N | vdhoeven2006; Edgar — |
| T. Hilberdink 2006 | Orders of growth of real functions (Bibliografie tid 1681) | Wachstumsordnungen | — | N | — |
| P. B. Miltersen u. a. 1999 | COCOON 1999 | halbexponentielle Schaltkreisgrössen | Komplexitätstheorie | S | miltersen1999 |
| W. Paulsen, S. Cowgill 2021 | Incompatible growth rates of the tetrations, Indian J. Math. 63(1), 127–141 | Kneser-Tetrationen verschiedener Basen liegen nicht in einem gemeinsamen Hardy-Körper; ihre gebrochenen Iterierten sind unverträglich | Hardy-Körper | S (zbMATH 1470.39046) | neu: bw.satz.paulsen_cowgill_hardy (siehe KARTE_bw) |
| P. Lévy 1927/1928; I. N. Baker 1955/1958 | gebrochene Iteration von eˣ über eˣ−1; keine ganze halbe Iterierte von eᶻ | Iteration | — | S (über L1) | levy1927, levy1928, baker1958 (lit_L1.yaml) |
| MathOverflow 4347 u. a. | „f(f(x)) = exp(x) and other functions just in the middle …“ | Frage nach halber exp und Wachstum | — | N | kn.satz.halbe_exp_wachstum (Forum-Quellen) |

## 4. Höhere Hyperoperationen

| Autor(en), Jahr | Werk, Ort | Hauptergebnis | Methode | P | Bei uns |
|---|---|---|---|---|---|
| A. A. Bennett 1915 | Note on an operation of the third grade, Ann. of Math. 17 | kommutative Hyperoperationen | — | N | bennett1915, br.satz.hyperop_varianten |
| W. Ackermann 1928 | Zum Hilbertschen Aufbau der reellen Zahlen, Math. Ann. 99 | nicht primitiv-rekursive Funktion φ(a,b,n) = Hyperoperationen | Rekursion | N | br.satz.ackermann (Quelle fehlt dort: Lücke an L3) |
| R. L. Goodstein 1947 | Transfinite ordinals in recursive number theory, J. Symbolic Logic 12 | Name „Tetration“, Folge der Hyperoperationen | — | S | goodstein1947 |
| J. Doner, A. Tarski 1969 | An extended arithmetic of ordinal numbers, Fund. Math. 65 | Hyperoperationen auf Ordinalzahlen | Ordinalarithmetik | N | — (L3) |
| D. E. Knuth 1976 | Coping with finiteness, Science 194 | Pfeilschreibweise | — | N | knuth1976 |
| K. A. Rubtsov, G. F. Romerio 1989/2004/2006 | Ackermann's function and new arithmetical operations (2004); Hyperoperations as a tool for science and engineering (ICM 2006, Poster); NKS-Fortschrittsbericht 2006/2007 | Zeration (Stufe 0), Hyperoperationen als Werkzeug; Fixpunkt der Tetration L_e ≈ −1.850354529 (nach MphLee, Forum tid 949) | Rekursion, Numerik | S (Forum tid 122, 949) | — (Lücke C; Zeration ggf. L3) |
| M. Müller 1995 | Reihenalgebra: What comes beyond exponentiation (mpmueller.net) | Fortsetzung der Hyperoperationen | — | N | — |
| H. Trappmann 2007 | Arborescent numbers: higher arithmetic operations and division trees | kommutative/Baum-Hyperoperationen | Algebra | N | br (L3) |
| D. Kouznetsov 2014 | Evaluation of holomorphic Ackermanns, Appl. Comput. Math. 3(6), 307–314, doi:10.11648/j.acm.20140306.14 | holomorphe Pentation zur Basis e („natürliche Pentation“, periodisch, Periode aus der Ableitung der Tetration am Fixpunkt ≈ −1.85), Algorithmen und Karten | reguläre Superfunktion der Tetration am reellen Fixpunkt | S (Abstract; Forum tid 949) | neu als Quelle für kn.konstr.pentation |
| D. Kouznetsov 2014/2020 | Суперфункции (2014, ISBN 978-3-659-56202-0); Superfunctions (LAP 2020, ISBN 978-620-2-67286-3) | Monografie: Superfunktionen, Abel-Funktionen, Tetration, natürliche Pentation | Übersicht | S (Verlagsangaben) | — |
| J. D. Fox 2007 | Forum tid 103 | Pentation zur Basis e (Matrix-slog), Werte auf ≈ 22 Stellen | Matrix-slog | O | kn.konstr.pentation, kn.beob.pentation_werte |
| S. Levenstein 2017 | Forum tid 1177 „pentation and hexation“ | Kette regulär/Kneser-artig: Pentation am linken reellen Fixpunkt der Tetration, Hexation Kneser-artig an den konjugierten Fixpunkten der Pentation (Kreis + θ, ≈ 30 Stellen), Heptation regulär; Schranken η₅ ≈ 1.738, η₆ ≈ 1.799, η₇ ≈ 1.839; hex(−3) = −2, Singularität nahe −3.86 | θ-Abbildung, Matrix | O | neu: kn.beob.hoehere_etas; Korrekturen an kn.beob.stufe5/6, kn.satz.hexation_struktur, kn.beob.unendliche_hexation |
| S. Levenstein 2020 | Forum tid 1263 „2^^^1.5 = ?“ | 2↑↑↑1.5 ≈ 2.7029, 2↑↑↑↑1.5 ≈ 2.6729 | fatou.gp/pentation.gp | O | forum_levenstein2020_pentation |
| H. Trappmann (nach Nixon 2017, tid 1177) | Forumsbeweis | jede analytische Hyperoperation der Stufe n erfüllt f(−k) = 1−k für 0 ≤ k ≤ n−1 | Induktion über Wertebereiche | Skizze (Beweis nicht gefunden) | S | kn.satz.hexation_struktur (Quelle ergänzen) |
| J. D. Nixon 2015/2021 | The bounded analytic hyper-operators, arXiv:1512.00831 (v2 zurückgezogen); A family of bounded and analytic hyper-operators, arXiv:2106.03935 | analytische, beschränkte Hyperoperationen aller Stufen für 1 < α < η | gebrochene Analysis (Mellin) | teils bewiesen | S | — (neu: Bezug in kn.verm.tpid10) |
| J. D. Nixon 2021 | Hyper-operations by unconventional means, arXiv:2101.03021 | Hyperoperationen über η-Basis-Iteration | — | S | — |
| A. Leonardis, G. d'Atri, F. Caldarola 2019; Caldarola u. a. 2019 | Beyond Knuth notation for unimaginable numbers; What are the unimaginable numbers (arXiv 1901.05372) | Notation sehr grosser Zahlen | Kombinatorik | S | — (L3) |
| J. D. Jaramillo 2020 | Hyperoperations in exponential fields, arXiv:2012.04093 | Tetration in exponentiellen Körpern (Modelltheorie) | Modelltheorie | S | — (L3) |
| J. P. Aguilera, A. Freund, M. Rathjen, A. Weiermann 2020 | Ackermann and Goodstein go functorial | Ackermann/Goodstein funktoriell | Beweistheorie | N | — (L3) |
| MphLee 2015 | Math.SE „On fractional ranks“ | gebrochene Stufen | — | N | — |
| Catullus 2022 | Forum tid 1532, 1594 | „schwache“ Hyperoperationen und 2sinh: eigene Etas (schwache Pentation ≈ 1.584 / 2.439) | — | S | neu in kn.frage.andere_etas |

## 5. TPID 10 — „Eulers → 4, Etas → 2“: alle gefundenen Zugänge

Bezeichnung hier: η₃ = η = e^{1/e}, e₃ = e; ηₙ, eₙ = Berührbasis und Berührstelle von x ↦ b↑ⁿ⁻²x mit y = x (η₄ für die Tetration usw.).

| Wer, wann | Wo | Zugang | Ergebnis | P |
|---|---|---|---|---|
| A. Robbins 11.04.2009 | Forum tid 269 „Infinite Pentation (and x-srt-x)“ | **Interpolation** der Selbst-Superwurzel srt_x(x) aus ganzen Stellen (Lagrange, mehrere hundert Punkte); keine Tetrations-Fortsetzung nötig | Maximum bei x ≈ 3.0885325, Wert ≈ 1.6353245; exakte Punkte 1.5596104695 (= ln2/W(ln2)) ↑↑↑∞ = 2, 1.6350784746 ↑↑↑∞ = 3 | O (Zusammenfassung) |
| nuninho1980 April 2009 | Forum tid 262 „tetration limit ??“, Beiträge 2, 11, 14 | **Matrix-slog** nach Robbins, 30 und 200 Dimensionen | 1.6353244967↑↑↑∞ ≈ 3.0885549441 (200 Dim.); Frage „Super-Euler?“ | O |
| bo198214 06.04.2009 | ebd., Beitrag 13 | Kennzeichnung durch slog_b(x) = x; **Methodenabhängigkeit** betont | — | O |
| mike3 03.07.2010 | Forum tid 468 | **Newton-Reihe** der Selbst-Tetrawurzel an a = 1 | Maximum bei ≈ 3.089, Wert ≈ 1.635; Pentation: Grenze zwischen beschränkt und explodierend | O |
| andydude 07.05.2010 | Forum tid 446 „Regular pentation?“ | Konvergenzintervall der unendlichen Pentation | obere Ecke (1.6353, 3.0885); **untere Grenze unbekannt** | O |
| S. Levenstein 31.10.2010 | tid 262, Beitrag 28 (pid 5327) | **Kneser** (kneser.gp), 32 Stellen | η₄ = 1.6353244967152763993453446183062, e₄ = 3.0885322718067176544821807826411; unterer Fixpunkt −1.64087257571659…, Multiplikator 4.80600575430175…; Taylorreihen | O |
| nuninho1980 31.10.2010 | tid 262, Beitrag 29 | Nachrechnung mit kneser.gp | 3.08853227181… bei der auf 23 Stellen genauen Basis (11 Stellen gleich) | O |
| dantheman163 31.10.2010 | tid 162, TPID 10 | **Vermutung** | Eulers → 4, Etas → 2 | O |
| tommy1729 01.12.2010 | tid 162, TPID 11 | **Vermutung (Ungleichung)** | ηₙ² < e_{n−1}; (ηₙ² − e_{n−1})/(η_{n−1}² − e_{n−2}) → 1 | O |
| Cherrina_Pixie, Levenstein, bo198214 Juni/Juli 2011 | Forum tid 668 „Base 'Enigma'“ | Pentation **an** η₄ selbst; parabolischer Fixpunkt e₄ ⇒ die reelle Pentation nimmt den unteren Fixpunkt; parabolische Fassung wie „cheta“ wäre anders | Pentation bei η₄: Fixpunkt ≈ 2.00, Periode ≈ 4.0024 i (beinahe ganzzahlig, unerklärt) | O |
| J. D. Nixon 09.07.2011 | Forum tid 675 „the infinite operator“ | **Heuristik** „unendlicher Operator“: 2 △_∞ 2 = 4 | Grenzwerte 2 und 4 als Folge davon, dass 2 ↑ⁿ 2 = 4 | O |
| D. Kouznetsov 01/2015 | Forum tid 949 | nicht-periodische Superfunktion („SuTet“) der Tetration zur Basis τ ≈ 1.63532 vorgeschlagen | — | S |
| S. Levenstein 03.09.2017 | Forum tid 1177, Beitrag 5 (pid 8759) | **höhere Glieder** mit Kneser-artiger Hexation und regulärer Pentation/Heptation | Schranken η₅ ≈ 1.738, η₆ ≈ 1.799, η₇ ≈ 1.839 (Funktionen auf 25–30 Stellen); Frage: beschränkt für alle b < 2? | O |
| J. D. Nixon 03./04.09.2017 | tid 1177, Beiträge 6, 8 | **Heuristik**: b > ηₙ ⇔ b↑ⁿx → ∞; daraus: Basen < 2 bleiben für grosse n beschränkt; Konstanten hängen von der gewählten Fortsetzung ab | Vermutung | O |
| S. Levenstein 04.09.2017 | tid 1177, Beitrag 7 | Frage „warum gehen die Eulers gegen 4?“; Basis 1.84: Heptation fast mit oberem Fixpunkt nahe 3.7, oct(45) ≈ 4.05 | — | O |
| Hyperops-Wiki, Stand 12.12.2022 | Seite „Euler number“ | drei Kennzeichnungen auf Stufe 4 (unendliche Pentation; Maximum der Selbst-Tetrawurzel; Fixpunktzahl 3 → 1) mit Zuschreibung (Andrew/Nuninho tid 269, Mike tid 468, Sheldon tid 262) | — | O |
| Catullus, Nixon, MphLee Juni 2022 | Forum tid 1518 „Natural properties of the tetra-Euler number“ | Ist e₄ „das nächste η“? Kennzeichnungen gesucht; Vorschlag allgemeiner Definition | Fragen | O |
| Catullus Juni/Juli 2022 | tid 1532, 1594, 1588 | Etas anderer Familien (schwache Pentation ≈ 1.584 / 2.439; 2sinh-Methode: Trappmann nennt e^{1/2} als „η“ von 2sinh); Ist e₄ rational? | Fragen | O/S |
| J. D. Nixon 2015/2021 | arXiv 1512.00831, 2106.03935 | beschränkte analytische Hyperoperationen aller Stufen für 1 < α < η | bewiesen (teils) | S |
| diese Arbeit 08.–09.10.2026 | exp.stufe5, exp.hexation_kneser | η₅, e₅, η₆, e₆ auf 30 Stellen; Abstandsgesetz 2 − ηₙ ≈ 2.8 n^{−1.47} | numerisch | O |
| OEIS (Suche 12.10.2026) | Ziffernfolgen 1,6,3,5,3,2,4,4,9,6,7 und 3,0,8,8,5,3,2,2,7,1,8 | kein Eintrag für η₄, e₄ (oder höhere Glieder) | — | O |
| MathOverflow, Math.SE (Websuche 12.10.2026) | — | keine Behandlung gefunden | — | S |

## 6. Die TPID-Liste (Forum tid 162, vier Seiten) und ihre Zuordnung

| TPID | Wer, wann | Inhalt | Stand laut Forum | Bei uns / Strang |
|---|---|---|---|---|
| 1 | bo198214, 17.05.2008 | Eigenwerte abgeschnittener Carleman-Matrizen von bˣ → (ln a)^k | offen | kn.verm.carleman_eigenwerte (L2) |
| 2 | andydude, 26.05.2008 | Eindeutigkeit des Exponentialfakultät; EF(0) | offen | ar (L3) |
| 3 | bo198214, 29.06.2008 | dasselbe für (x+s)^p − s, Grenzwerte p^k | offen; in voller Allgemeinheit numerisch zweifelhaft | kn.verm.tpid3, kn.verm.tpid3_bedingung (exp.tpid3) |
| 4 | bo198214, 08.10.2008 | keine ganze Lösung von f(z+1) = b f(z) beschränkt im Streifen | offen; tommy1729 2022 Skizze | fg.satz.tpid4 (L1) |
| 5 | andydude, 01.05.2009 | √2, f(−f(x)) = −x ⇒ nirgends stetig | offen | kl.frage.tpid5 (L1) |
| 6 | andydude, 07.10.2009 | Lösung von ⁿx = n → η | **bewiesen** (Forum) | fp.satz.tpid6 (L1) |
| 7 | andydude, 23.10.2009 | Schranken für (^y x)(^{−y} x) | offen | gr.verm.tpid7 (L3), kl (L1) |
| 8 | bo198214, 25.04.2010 | elementare F mit F(1+F⁻¹(x)) Polynom ohne reelle Fixpunkte | offen | it (L1) |
| 9 | bo198214, 20.07.2010 | Folge aₙ → (b−1)/ln b | **bewiesen** (Trappmann, arXiv:1008.1409) | kn.satz.intuitiver_log (L2) |
| 10 | dantheman163, 31.10.2010 | Eulers → 4, Etas → 2 | offen | kn.verm.tpid10 (L2) |
| 11 | tommy1729, 01.12.2010 | ηₙ² < e_{n−1}; Quotient → 1 | offen | bisher nur Satz in kn.verm.tpid10 — neu kn.verm.tpid11 (L2) |
| 12 | bo198214, 31.05.2011 | Interpolationspolynome der Superwurzel-Punkte | offen | fp.frage.tpid12_13 (L1) |
| 13 | bo198214, 31.05.2011 | Interpolation von n^{1/n} → x^{1/x}? | offen | fp.frage.tpid12_13 (L1) |
| 14 | tommy1729, 01.06.2011 | Robbins' slog: Entwicklungspunkt x₀ durch 0 ersetzen gibt denselben slog | offen | kn.frage.tpid14 |
| 15 | tommy1729, 27.05.2014 | Lösung von x↑↑x = −1 ist nicht negativ reell | Bemerkung | — (L3, gr/ar) |
| 16 | tommy1729, 07.06.2014 | ganze reelle Funktionen mit konjugiertem primärem Fixpunktpaar, Positivität der Ableitungen der Iterierten | offen | — Lücke C (L2 oder L1) |
| 17 | tommy1729, 28.03.2015 | Fehlerschranke O(n ln(n+1)) für S9 der „fake function theory“ (tid 863) | offen | — Lücke C (L2, halbe exp) |
| 18 (a) | andydude, 25.12.2015 | Reihe für die dritte Superwurzel (tid 1033) | offen | gr.satz.superwurzel_dritte (L3) |
| 18 (b) | Gottfried (G. Helms), 04.07.2022 | „Polynom-Tetration“ aus abgeschnittenen Carleman-Matrizen → Kneser? | offen | kn.verm.carleman_kneser (f3), Doppelnummer im Forum |
| — | Catullus 12.07.2022; Leo.W 10.08.2022; Catullus 01.11. und 22.12.2022 | unnummerierte Vermutungen (exp einzig total monoton in der Basis; Fixpunkte/Multiplikatoren; van-der-Corput-Grenzwerte) | offen | kl.verm.vollmonoton (L1) u. a. |

## 7. Lückenliste kn (Priorität A = zuerst)

- **A (erledigt 12.10.)** Eigene „vermutlich neu“-Einstufungen bei Hexation/höheren Stufen korrigieren: η₅, η₆ (Levenstein 2017 auf 3–4 Stellen), unendliche Hexation bis η₅ (ebd.), H(−3) = −2 und Singularität (ebd.; Trappmann nach Nixon), Kneser-artige Hexation über Kreis und θ (ebd.). → erledigt in diesem Durchgang (BRETT).
- **A (erledigt 12.10.)** TPID 10 mit allen Zugängen (Abschnitt 5) als Knoten (Strang kn.etas, 9 Aussagen); TPID 11 als eigene Vermutung.
- **A (wartet auf den Autor)** Volltexte Paulsen 2026 (z↑↑a als Funktion der Basis, betrifft bw.ann.basis_diffbar) und Paulsen–Cowgill 2021 — nicht frei zugänglich; bis dahin sekundaer.
- **A (zurückgestellt, Entscheid des Autors)** η₇, e₇ genau rechnen (Levenstein 2017: η₇ ≈ 1.839; Abstandsgesetz 1.8399; e₇: Hochrechnung ≈ 3.54 gegen Levensteins „nahe 3.7“).
  Aufwand: neue Stufe — reguläre Heptation am linken reellen Fixpunkt der Kneser-artigen Hexation (Schröder-Reihe der Hexation dort, Hexation auf einem Kreis
  um den Berührpunkt). Code etwa ein Arbeitstag; Rechnung etwa 10–15 Basen in der Sekantensuche, je Basis ein knhexinit-Lauf plus einige hundert Hexationswerte —
  geschätzt 1–4 h Laptop, Kneser-Motor exklusiv (zum Vergleich: exp.hexation_kneser 141 s mit Zwischenspeicher für eine Basis).
- **B (erledigt 12.10.)** Kouznetsov 2014 an kn.konstr.pentation; Kouznetsov 2010 als Quelle an kn.konstr.cauchy.
- **B (erledigt 12.10.)** Szekeres 1998 → kn.verm.szekeres_wachstum.
- **B (erledigt 12.10.)** TPID 3 (mit Befund: gilt nur, wenn 1 − s im Konvergenzkreis liegt) und TPID 14 als Knoten.
- **B** Paulsen 2019 im Volltext: zeigt es die Holomorphie in b (Paper V, (K3))? Bis dahin als „nicht belegt“ markiert.
- **C** Bromer 1987 (Superexponentiation, Math. Mag. 60, 169–174) und Kouznetsov–Trappmann 2010 (Superfunctions and square root of factorial, Mosc. Univ. Phys. Bull. 65, 6–12) — nur kartiert.
- **C (teilweise)** Clenshaw–Lozier–Olver–Turner 1986 als Quelle an kn.konstr.linear (Definition nicht gelesen); Tetration FAQ als Sekundärquelle der Verfahren offen.
- **C** Hyperoperations-Literatur (Rubtsov–Romerio, Müller, Doner–Tarski, Ackermann 1928) — Grundlagen gehören zu br (L3); hier nur Kartierung.
- **C** TPID 16, 17 (tommy1729) einordnen.

## 8. Stand und Restlücken (12.10.2026, nach Block 3)

**Stand:** Kartiert sind Kneser und seine Eindeutigkeit (Kouznetsov 2009 → Trappmann–Kouznetsov 2011 → Paulsen–Cowgill 2017), die anderen Verfahren
(Walker, Robbins/Carleman, Szekeres, Paulsen 2016, Kontinuumssumme, 2sinh, β-Methode), die halbe exp, die höheren Hyperoperationen (Kouznetsov 2014, Levenstein 2017/2020),
TPID 10 mit allen Zugängen (Strang kn.etas) und die ganze TPID-Liste: TPID 1, 3, 9, 10, 11, 14, 16, 17, 18 (b) haben Knoten im Gebiet kn. Robbins 2005: Verfahren nach Forum tid 3
belegt (Arbeit selbst nicht gelesen; Trappmann: Ansatz schon bei Walker um 1990). KANDIDATEN-Zeile Paulsen 2019 = `paulsen2019` (Adv. Comput. Math. 45, 243–267, über zbMATH bestätigt).

**Restlücken:**
- Volltexte hinter Bezahlschranke (bleiben sekundaer): Paulsen–Cowgill 2017 und 2021, Paulsen 2019 (Holomorphie in b, Paper V (K3) — „nicht belegt“), Paulsen 2026, Kneser 1950 (nur über Forum tid 213), Kouznetsov 2009/2010, Walker 1991.
- η₇, e₇ genau (zurückgestellt; Aufwand §7).
- Nicht erreichbar: Math.SE (Abruf gesperrt); Robbins' Arbeit (tetration.itgo.com) nicht gelesen.
- Nur kartiert, ohne Knoten: Bromer 1987 (Math. Mag. 60, 169–174, doi:10.1080/0025570X.1987.11977296), Wassell 2000 (Math. Mag. 73, 111–119; Fixpunkte → L1),
  Kouznetsov–Trappmann 2010 (Fakultät), Clenshaw u. a. 1986 (Definition nicht gelesen), Tetration FAQ, Ueda 2021, Nesargi–Roudenko 2025, Aldrovandi 2014;
  Hyperoperations-Grundlagen (Rubtsov–Romerio, Müller, Doner–Tarski, Ackermann 1928) gehören zu br (L3).

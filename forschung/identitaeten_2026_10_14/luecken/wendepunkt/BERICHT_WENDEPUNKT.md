# Wendepunkt von Knesers sexp: Zertifikat mit Intervallarithmetik (14.10.2026)

## Kurz gesagt

- **Ergebnis:** Knesers sexp hat zur Basis e genau einen Wendepunkt auf (−2, ∞). Er liegt bei
  **x_w = −0.521764792623322155340028093 ± 4·10⁻²⁸**. Zur Basis 2 gibt es ebenfalls genau einen, bei
  **x_w = −0.030244542753639263604287235 ± 3·10⁻²⁸**. Beide Nullstellen von sexp″ sind einfach.
- **Was streng gerechnet ist:** Jeder Schritt von den Motorwerten bis zur Aussage läuft in Ballarithmetik (python-flint arb/acb).
  Rundungsfehler, Abbruchfehler der Taylor-Reihen, Aliasing der diskreten Fourier-Transformation und die Kettenregel über die
  Turmgleichung sind alle als garantierte Schranken eingerechnet.
- **Die einzige ungeprüfte Annahme H:** Die Motorwerte (fatou.gp) von sexp an den 5 × 256 Knoten der Kreise |z − c| = 0.2 sind
  auf **10⁻⁵** genau (H1), und es gilt |sexp| ≤ **1000** auf den Kreisscheiben |z − c| ≤ 0.3 (H2), mit c ∈ {−0.4, −0.2, 0, 0.2, 0.4}.
  Das ist eine sehr schwache Annahme. Gemessen stimmen die Werte bei verschiedener Motorgenauigkeit auf 10⁻⁵¹ (Basis e) bzw.
  10⁻⁷⁰ (Basis 2) überein, und max |sexp| auf den Kreisen ist ≈ 2. Mit der strengeren Annahme 10⁻³⁰ / C = 10 ergibt sich die
  oben genannte Genauigkeit von x_w. Unter 10⁻⁵ reicht die Genauigkeit nur noch für die Aussage selbst, x_w dann auf ±4·10⁻³.
  Bei 10⁻³ schlägt das Zertifikat fehl (Gegenprobe: die Prüfung ist nicht leer).
- **Ist es wirklich Kneser?** Das lässt sich nicht streng zeigen, ohne den Motor selbst zu zertifizieren. Belegt ist aber
  Folgendes: Die Motorwerte erfüllen die Turmgleichung (10⁻⁹⁰). Sie stimmen auf 3·10⁻⁵⁹ mit sexp = F(z + θ(z)) überein, wobei die
  θ-Reihe nur nicht-negative Fourier-Moden hat. Das ist Knesers Kennzeichen nach Paulsen–Cowgill. Wer H1 unabhängig nachprüfen will,
  braucht nur 5 Stellen einer beliebigen Kneser-Implementierung.
- **Anders als geplant:** Für die reelle Aussage braucht man die Windungszahlen im Komplexen gar nicht. Es genügt, das Vorzeichen
  von sexp″ auf einem reellen Stück zu zertifizieren: [−1.6, 0.6] für e und [−1.1, 1.0] für 2. Den Rest liefern zwei elementare
  Lemmata über die Turmgleichung. Damit ist auch kn.satz.ableitung_nullstellenfrei hier nicht nötig. Die *komplexe* Aussage W
  („sexp″ hat in [−1.8, 1.4] × [−4, 4] nur die eine Nullstelle“) ist **nicht** zertifiziert und bleibt rechnergestützt.

## 1. Der Beweis (b > η, T = Knesers sexp_b, λ = ln b)

Aus T(x+1) = e^{λT(x)} folgen exakt:

- (L1) T′(x+1) = λ T(x+1) T′(x), und T″(x+1) = λ T(x+1) · (λ T′(x)² + T″(x)).
- (L2, nach rechts = gr.kor.konvex_nach_oben) Ist T″ ≥ 0 und T′ > 0 auf einer ganzen Periode [x₀, x₀+1], so ist T″ > 0 auf
  [x₀+1, ∞). Das folgt aus (L1) mit Induktion, denn T > 0 rechts von −1.
- (L3, nach links) Es gilt T″(x) = T″(x+1)/(λ T(x+1)) − λ T′(x)². Ist T″(x+1) < 0 und T(x+1) > 0, so ist T″(x) < 0.

**Zertifiziert, Basis e** (`w_reell.py e`):
- T′ > 0 auf [−1.6, 0.6].
- T″ < 0 auf [−1.6, −0.528].
- T‴ > 0 auf J = [−0.528, −0.516].
- T″ > 0 auf [−0.516, 0.6].

Daraus folgt: T″ ist auf J streng wachsend, an den Rändern von J negativ bzw. positiv, hat also in J genau eine Nullstelle, und die
ist einfach. Intervall-Newton schliesst sie ein auf x_w ∈ [−0.521764792623322155340028093 ± 3.96·10⁻²⁸].
- Für [0.6, ∞): T″ > 0 und T′ > 0 auf der ganzen Periode [−0.4, 0.6], also gilt (L2).
- Für (−2, −1.6): Dann liegt x+1 in (−1, −0.6) ⊂ [−1.6, −0.528]. Dort ist T″ < 0. Ausserdem ist T(x+1) > 0, weil T(−1) = 0 und T′ > 0.
  Mit (L3) folgt T″(x) < 0.

Also hat T″ auf (−2, ∞) genau eine Nullstelle, und an ihr wechselt das Vorzeichen. ∎ (unter H)

**Basis 2** (`w_reell.py 2`):
- T′ > 0 auf [−1.1, 1.0].
- T″ < 0 auf [−1.1, −0.032].
- T‴ > 0 auf [−0.032, −0.030].
- T″ > 0 auf [−0.030, 1.0].

Die Nullstelle liegt bei x_w ∈ [−0.030244542753639263604287235 ± 2.34·10⁻²⁸]. Rechts von 1 gilt (L2) mit der Periode [0, 1]. Links,
auf (−2, −1.1), gilt (L3), denn x+1 liegt in (−1, −0.1). ∎ (unter H)

Werte am Wendepunkt (Basis e): sexp(x_w) = 0.477743094766666235175664095…, sexp′(x_w) = 0.956521797528063948666374895…,
sexp‴(x_w) = 0.98882753092277431206….
Basis 2: sexp = 0.973106860339015173919356826…, sexp′ = 0.8891021633099558803459979568…, sexp‴ = 0.57632755562292080931….

## 2. Vom Motorwert zur Ballschranke (alles streng, `w_reell.py`)

Für jedes Zentrum c (r = 0.2, M = 256, R = 0.3, K = 128 Glieder):
1. **Knoten.** Der Motor bekam z̃_j, die Dezimaldarstellung von z_j = c + r e^{2πij/M} mit 140 bzw. 90 Stellen. Es gilt
   |sexp(z̃_j) − sexp(z_j)| ≤ |z̃_j − z_j| · C/(R − r − 10⁻³) (Cauchy). Jeder Wert wird zum Ball s_j mit Radius
   ε + diese Verschiebung.
2. **Diskrete Fourier-Transformation in acb.** Es gilt (1/M) Σ_j sexp(z_j) e^{−2πijk/M} = Σ_{m≥0} a_{k+mM} r^{k+mM}.
   Die Aliasglieder (m ≥ 1) sind höchstens C R^{−k} q/(1−q) mit q = (r/R)^M ≈ 10⁻⁴⁵. Sie werden als Radius addiert.
   a_k ist reell, weil sexp auf ℝ reell ist. Der Imaginärteil-Ball enthält jedes Mal 0; das ist eine Probe.
3. **Reihenrest** für die j-te Ableitung (j ≤ 3) auf |x − c| ≤ ρ = 0.11: C Σ_{k≥K} k³ ρ^{k−3} R^{−k} ≤ 4·10⁻⁴⁶ (bei C = 10).
4. **Auswertung.** Das Polynom mit Ballkoeffizienten wird auf Teilintervallen der Breite 1/500 ausgewertet, mit Halbierung bei
   Bedarf (es waren nie welche nötig). Ausserhalb [−½, ½] geht es über T(x) = log T(x+1)/λ bzw. T(x) = e^{λT(x−1)} und die
   Kettenregel bis zur dritten Ableitung, ebenfalls in Bällen.
5. Die Vorzeichenlogik und Intervall-Newton (N(I) = m − T″(m)/T‴(I), geschnitten mit I) laufen auf exakten Brüchen bzw. Bällen.

Laufzeit je Basis etwa 2 s, Arbeitsgenauigkeit 640 Bit.

## 3. Belege für die Annahme H (Gleitkomma, `w_pruef.py`, Protokoll `log_pruef.txt`, `pruefung.json`)

| Prüfung | Ergebnis |
|---|---|
| P1 gleiche Knoten, Motor dps 130 gegen 100 (Basis e) | ≤ 5.5·10⁻¹¹⁰ |
| P1 dps 130 gegen 70 (Basis e) | ≤ 2.1·10⁻⁵¹ |
| P1 dps 80 gegen 60 (Basis 2; dps 60 hier neu gerechnet, Speicher `cache/2_60.json`) | ≤ 5·10⁻⁷⁰ |
| P2 Turmgleichung: exp(Motor auf dem Kreis um c−1, dps 80) gegen Motor auf dem Kreis um c (dps 130), Knoten für Knoten, 5 Paare | ≤ 5.4·10⁻⁹⁰, keine Ausreisser |
| P3 Kneser-Kennzeichen: F(z + θ(z)) mit θ aus Etappe 2 (nur Moden n ≥ 0), gegen den Motor an allen Knoten mit Im z ≥ 0.1 (je 43) | e: ≤ 3.6·10⁻⁵⁹, 2: ≤ 4.7·10⁻⁵⁹ (Grenze ist die Rechengenauigkeit dps 60) |
| P4 max \|sexp\| auf \|z − c\| = 0.3 (θ-Reihe bzw. Taylor-Daten) | e: 2.00, 2: 1.66 (Annahme H2: ≤ 10 bzw. 1000) |

**Empfindlichkeit** (`python w_reell.py e EPS C`):

| ε (H1) | C (H2) | Basis e | Basis 2 |
|---|---|---|---|
| 10⁻³⁰ | 10 | zertifiziert, x_w ± 4·10⁻²⁸ | zertifiziert, ± 2·10⁻²⁸ |
| 10⁻⁸ | 1000 | zertifiziert, ± 7·10⁻⁶ | zertifiziert, ± 6·10⁻⁶ |
| 10⁻⁵ | 1000 | zertifiziert, ± 4·10⁻³ | zertifiziert, ± 2·10⁻³ |
| 10⁻³ | 1000 | schlägt fehl (T‴ > 0 auf J nicht zeigbar) | – |

Nebenbefund: In `etappe3/daten/taylor_links_80_256.json` sieht a₀ beim Zentrum −1.0 wie „5.406…“ aus. Vollständig gelesen heisst
der Eintrag 5.4·10⁻⁸¹, ist also in Ordnung (sexp(−1) = 0). Die Knotenwerte dort bestehen P2.

## 4. Was nicht zertifiziert ist

1. **Die Identität der Daten mit Kneser.** H1 setzt voraus, dass der Motor Knesers sexp rechnet. P3 zeigt die Kneser-Struktur
   numerisch. Ein strenger Nachweis bräuchte eine zertifizierte θ-Reihe: eine Schranke |θ − θ̃| ≤ ε auf einer Geraden Im z = y₀,
   die dann nach dem Maximumprinzip für alle Im z ≥ y₀ gilt, und dazu Realität am Rand. Das ist nicht gemacht.
2. **Die komplexe Nullstellenzählung W** auf [−1.8, 1.4] × [−4, 4] und |φ − L| ≪ |L| oberhalb von Im = 4. Für die reelle
   Wendepunkt-Aussage ist sie überflüssig. Als eigene Aussage (keine komplexen Nullstellen von sexp″ in der Nähe) bleibt sie
   rechnergestützt. Für eine Zertifizierung bräuchte man Bälle für θ, θ′, θ″ aus einer θ-Hypothese wie in Punkt 1, Ψ mit Restschranke
   und Argumentzuwächse mit dem Kriterium |f(z_{k+1}) − f(z_k)| < |f(z_k)|.
3. **Nur zwei Basen.** Die Familie in b (gr.verm.ein_wendepunkt) ist offen. Für jede weitere Basis reichen die Motorwerte auf fünf Kreisen und
   ein Lauf von `w_reell.py`. Wo der Wendepunkt weit rechts liegt (b nahe η), braucht man ein längeres Prüfstück und mehr Zentren.

**Vorschlag für den Graphen** (nicht eingetragen, `daten/` unverändert):
- gr.beob.ein_wendepunkt_e: Status von „numerisch“ auf „zertifiziert unter Annahme H (Motorwerte auf 10⁻⁵, |sexp| ≤ 1000)“ setzen,
  den reellen Weg L1–L3 statt der Windungszahlen nennen und x_w mit 27 Stellen angeben.
- Basis 2 als eigene Beobachtung eintragen, x_w = −0.0302445427536392636….

## 5. Dateien und Aufrufe

Interpreter `C:/Users/janis/miniconda3/python.exe`, `PYTHONIOENCODING=utf-8`, Ordner `luecken/wendepunkt/`.

| Datei | Inhalt |
|---|---|
| `w_daten.py` | liest die Motorwerte der Kreise aus `etappe3/cache/` (nur lesen) |
| `w_reell.py e` / `w_reell.py 2` (optional `EPS C`) | das Zertifikat, je ~2 s → `zertifikat_<b>_eps<EPS>_C<C>.json` |
| `w_pruef.py` | P1–P4 (~20 s; Basis 2 dps 60 einmal über den Motor, Speicher `cache/`) → `pruefung.json`, `log_pruef.txt` |

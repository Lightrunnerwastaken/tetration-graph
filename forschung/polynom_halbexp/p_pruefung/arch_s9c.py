"""Hilfsskript (Strang P): ARCHITEKTUR §9 auf Stand Runde 3 bringen. Aufruf: python arch_s9c.py NN ZEIT"""
import sys
from pathlib import Path

NN = sys.argv[1]
ZEIT = sys.argv[2]
TMAX = sys.argv[3] if len(sys.argv) > 3 else "≈ 35 s (N ≈ 190)"
p = Path(__file__).resolve().parent.parent / "ARCHITEKTUR.md"
s = p.read_text(encoding="utf-8")


def ersetze(alt, neu, n=1):
    global s
    assert s.count(alt) >= 1, alt[:80]
    s = s.replace(alt, neu, n)


ersetze("Stand: 11.10.2026, 00:42 (Zeiten aus `date`). Einzelnachweise: `p_pruefung/PRUEFLISTE.md`, Skripte p01–p17.",
        f"Stand: 11.10.2026, 00:42; Runde 3 eingearbeitet {ZEIT} (Zeiten aus `date`). Einzelnachweise: `p_pruefung/PRUEFLISTE.md`, Skripte p01–p24.")

ersetze("| (A) | **Lücke (Kopplung), ausgearbeitet in §9.7** | R_N ≈ 240·C·N^{κ−4.81}; (A) = (H) ab N₀ + Zweig-Zertifikat bis N₀ (p18: rigoros bis N ≥ 45) + zertifizierte f_K auf D̄ auf ~1e-6 |",
        f"| (A) | **Lücke (Kopplung), ausgearbeitet in §9.7** | R_N ≈ 240·C·N^{{κ−4.81}}; (A) = (H) ab N₀ + Zweig-Zertifikat bis N₀ (**rigoros bis N = {NN}**: p18 bis 45, p24 kohärent ab 45) + zertifizierte f_K auf D̄ auf ≈ R_{{N₀}} (**fehlt; Weg §9.8**) |")

ersetze("""(Kreuzung bei N ~ 1e14). In W̄ konvergiert dasselbe Profil (0.22), aber der A-Bogen liegt in f_K(D) (Eingaben ~(ρ/|h|)^N, p02). K* nimmt
die Ecke von W̄ und den Rest von K'; ∂K* ∩ f_K(D) liegt dann in A_r := D(L,2r₁) ∪ D(L̄,2r₁) (A-Bogen-Stück nahe L plus Verbindungsbogen zu
f_K(Γout)).""",
        """(Kreuzung bei N ~ 1e14). In W̄ konvergiert dasselbe Profil (0.22), aber der A-Bogen liegt in f_K(D) (Eingaben ~(ρ/|h|)^N, p02). K* nimmt
die Ecke von W̄ und den Rest von K'; ∂K* ∩ f_K(D) liegt dann in A_r := D(L,2r₁) ∪ D(L̄,2r₁) (A-Bogen-Stück nahe L plus Verbindungsbogen zu
f_K(Γout)).
**KORREKTUR (P, Runde 3, nach W 01:06/§14):** Der Exponent 1.26–1.29 ist zurückgezogen (W1 zurückgezogen); das massgebliche Γin-Profil (Eckteil
der kopffreien Pol-Antwort, kubisch UND Kneser) hat **p = 1.5–1.6**. Mit Ws neuem Profil (λ = 0.0001/0.0003/0.0009/0.0016/0.0024 bei
0.5°/1°/2°/3°/4°) und konservativ λ = 0 auf der anderen Seite gibt p17b den Ecken-Grenzwert **K': 0.191 (p = 1.5) / 0.155 (p = 1.6)**, W̄: 0.089/0.080
— beide konvergent (p·α' = 1.18 > 1). **K' ist also NICHT grenzwertig**; der Hybrid K* und (E_A) werden für die Ecke nicht mehr gebraucht.
Der reine K'-Weg (P-II für den Rest, P-III mit p ≈ 1.5 für den Eckteil, bzw. Variante X1) ist damit gleichwertig zum gemischten Weg; X 01:08
misst zusätzlich Abklingen des kanonischen Eckteils auf der f_K(Γout)-Seite (Ecken-Grenzwert −30).""")

ersetze("""Achtung: (Ẽ) in K' hat dieselbe Grenzwertigkeit an der Ecke wie oben (Γin-Exponent ≈ 1/α'); sie verlangt zusätzlich Abfall auf der
f_K(Γout)-Seite (X: gemessen ≈ e^{−0.016N}), oder man nimmt für c_N Poisson–Jensen in Ĉ∖W̄ (dann braucht man c_N auf dem A-Bogen, Typ (V3)).""",
        """~~Achtung: (Ẽ) in K' grenzwertig~~ — überholt (Runde 3): mit p ≈ 1.5 ist die K'-Ecke nicht grenzwertig (p17b), X 01:08 misst Θ ≈ 0.13 für den
kanonischen Eckteil (N = 200…600, kleine d) bzw. Θ^c < 0 für d = 2N (N ≤ 300, X 01:22 = Satz X3). Für d = 2N zeigt P (PRUEFLISTE, X 01:22) den Schwanz
‖T_{>3N}J_∞(f_K)δ‖ ≤ CN^{−3.95}‖δ‖ (Darboux an L, L̄; unter der Bedingung (T3): max_{|x|=R}|f_K| < (R/ρ)³ für ein R < |L ± 2πi|).""")

ersetze("Satz P-III mit α = 0.883 und p ≈ 1.27 ⇒ (E_K^{∂*})",
        "Satz P-III mit α = 0.883 und p ≈ 1.5 (W §14; früher 1.27) ⇒ (E_K^{∂*})")

ersetze("Kantenschwanz ‖Q_NY_NP_{cN}‖ ≤ CN^b [W, nach O §16] · (Fam) · (A) [P, Zertifikat in Arbeit].",
        f"Kantenschwanz ‖Q_NY_NP_{{cN}}‖ ≤ CN^b [W, nach O §16] · **(Fam) = OFFENE AUFGABE F (9.6)** · (A) [P: Zweig rigoros bis N = {NN}; es fehlt das Kneser-Zertifikat §9.8].")

zeilen_95 = f"""| W 01:06/01:29 (§14) | Exponent 1.27 zurückgezogen, Profil p = 1.5–1.6; (E_L)-Weg trägt in W̄ (Θ_∞ 0.22/0.10; Kneser 0.215/0.103), (E_A) gestützt | **übernommen; Folgerung nachgerechnet** (p17b: K'-Ecke 0.19/0.155, nicht grenzwertig — §9.3 korrigiert) |
| X 01:08 | kanonischer Eckteil in K': Θ ≈ 0.13, Ecken-Grenzwert −30 | plausibel, konsistent mit p17b (nicht nachgerechnet) |
| X 01:22 | Satz X3 (d = 2N: (a″)_∞ + einteilige Randbedingung Θ^c), Kneser Θ^c < 0 (N ≤ 300) | Reduktion wie X1 (verifiziert); Schwanz ‖T_{{>3N}}J_∞δ‖ ≤ CN^{{−3.95}}‖δ‖ von P skizziert (Darboux, Bedingung (T3)) |
| O 01:11 | O7g (Kommutator), O7h (Spuren, Carleson), O7i (Normstetigkeit) ⇒ O7 bis auf Standardsätze; O7j (T invertierbar für Kneser, qualitativ); o22: σ_min(T_{{3N}}J_∞|P_N) > 0.7315, N ≤ 100 | O7g/O7h/O7i **verifiziert** bis auf CMM/David, Calderón, Pommerenke, Carleson, Kato (Standard); O7j **verifiziert|(G2),(S),P1–P4,D4** (D4 nur Skizze); o22 **bedingt** (Kneser-Koeffizienten nicht zertifiziert, §9.8) |
| M 01:04/§8 | keine Gårding-Ungleichung; (a″)_d mit d ≈ r(N) ≡ (a′)_r; Lemma M3 (σ_min ≤ min_Γin|m| + C/N für alle d — 0.690525 scharf); Lemma M2 | M3 **verifiziert bis auf F = 0.8751 < 1 (Gitter)**; M2 **verifiziert** |
| P Runde 3 | **Kohärentes Zweig-Zertifikat** (p24): jeder Übergang N → N+1, N = 45 … {int(NN) - 1}, mit EINEM s-Segment | **rigoros** (arb), §9.7 (3b) |
| P Runde 3 | Kneser-Werte mit Schranken: in gp-tetration NICHT vorhanden (nur Vorzertifikate); eigener Weg Satz P-V | **Befund + Bauplan**, §9.8 |
"""
ersetze("\n### 9.6 Abhängigkeiten (Stand 00:42)", "\n" + zeilen_95.rstrip("\n") + "\n\n### 9.6 Abhängigkeiten (Stand 00:42; Runde 3 ergänzt)")
# Die neuen Zeilen gehören in die Tabelle von 9.5: direkt vor der Leerzeile einfügen
s = s.replace("\n\n" + zeilen_95.rstrip("\n") + "\n\n### 9.6", "\n" + zeilen_95.rstrip("\n") + "\n\n### 9.6")

ersetze("""(H₀*) bzw. (H₀^mix) ⇒ (H) [Familien-Stabilität, offen] · alternativ einteilig: R_N(L) polynomial + Θ*_N ≤ 1−ε ⇒ (H₀*) [P-I(ii)].""",
        f"""(H₀*) bzw. (H₀^mix) ⇒ (H) [Familien-Stabilität, offen] · alternativ einteilig: R_N(L) polynomial + Θ*_N ≤ 1−ε ⇒ (H₀*) [P-I(ii)].

**OFFENE AUFGABE F — Familien-Stabilität (H₀) ⇒ (H)** (sichtbar geführt; die Leitung vergibt sie).
*Gegeben:* (H₀) in einer der Formen 9.3 (K', gemischt, X1/X3) für die feste Linearisierung J_∞ = J_∞(f_K).
*Zu zeigen:* (H-i)–(H-iii) aus §5.1 für ALLE Fast-Lösungen P ∈ 𝒮_N ∩ B_N, also für die Sekantenoperatoren 𝒥_P v = v∘P + M_P v mit
P − f = λ(Π^P_N f − f) + Kopfterm, sup_{{D̄}}|P − f| ≤ R_N = O(N^{{κ−4.69}}): dieselben Pol-Antwort-Schranken (bis auf einen Faktor) und die Kopf-Stabilität.
*Die Schwierigkeit:* v∘P − v∘f braucht v auf P(D̄) ≈ f_K(D̄) ⊄ D̄ — für Polynome vom Grad N nur mit e^{{N g}}-Verlust kontrollierbar; tragfähig sind nur
**kohärente** Störungen (Π^P f − f klein auf ganz K', nicht nur auf D̄). Eine Kugel aller Polynome um f in sup_D̄ ist zu gross (p22: die Kopfantworten
J_N⁻¹x^k sind auf P(D̄) exponentiell gross, ~1.3^N).
*Stützen:* p09 (λ-Pfad N = 60: Pol-Antworten unverändert, (H-iii)-Quotient 1.0000); p24 (längs der echten Übergänge ändert sich J_N relativ nur um
θ ≤ 6e-8, N ≤ {NN}). *Bausteine:* Störungsformel E^P_N − E_N über 𝒥_P − J_∞ = (Komposition) + (M_P − M_f); Os c_∞ > 0 ist offen in P (Kato).
*Entlastung durch p24:* (Fam) wird nur noch für N ≥ N₀ gebraucht, N₀ ≤ {NN} frei wählbar.""")

ersetze("""**Fazit:** (A) ist kein unabhängiger Test, sondern: [(H) mit expliziten Konstanten ab N₀] + [Zertifikat bis N₀] + [zertifizierte f_K-Werte auf D̄ auf ~1e-6].""",
        f"""**(3b) Kohärentes Zweig-Zertifikat (p24, Runde 3, rigoros, arb) — das 1.62^N-Wachstum ist beseitigt.** Statt einer Koordinaten-Box um den GANZEN
Pfad eine Röhre Z(s) = z̃(s) + B: z̃ = Interpolationspolynom (Grad 6) durch Newton-Punkte des Pfades an 7 Tschebyschow-Knoten in [0,1], B eine winzige
Box (Radien r·ρ^{{−k}}, r ≈ 1e-55…1e-60). Krawczyk für alle s ∈ [0,1] und c ∈ C zugleich:
K − z̃(s) = −Y G_{{s,c}}(z̃(s)) + (I − Y J(z̃(s)))B − Y(J(Z) − J(z̃(s)))B ⊂ int B. Dabei werden s ↦ G(z̃(s)) und s ↦ J(z̃(s)) (Polynome in s) NUR an den
Knoten in Punktarithmetik ausgewertet (kein Wrapping) und über die Tschebyschow-Restformel (|ω| ≤ 2^{{1−2n}}, Hermite–Genocchi) auf [0,1] fortgesetzt;
die n-te Ableitung wird über Majoranten-Reihen auf |s| ≤ R (R ≈ 1/(M‖â₁‖) ~ 1e8…1e11) mit Cauchy beschränkt. Das Wrapping 1.62^M steht nur noch
vor R^{{−7}} bzw. vor r (Präzision 2M + 320 Bit genügt; Knoten auf 2^{{−prec}} gerundet — |ω| ≤ 2^{{1−2n}} und Λ ≤ (2/π)ln n + 1 gelten
dann bis auf relative 1e-100, gedeckt durch die Reserve 0.5; Kopfkoeffizient c als Intervall C über die Box von p_N). Ergebnis: **alle Übergänge N → N+1 für N = 45 … {int(NN) - 1} mit EINEM s-Segment**, Reserve 0.5,
linearer Term θ = sup_s‖I − YJ(z̃(s))‖_w ≤ 6e-8 (fällt mit N), Wrapping-Term ≤ 3e-8 (relativ; wächst wie ~1.62^N·r — für N ≳ 230 braucht es
kleineres r, d. h. mehr Präzision und Knoten); Rechenzeit 1 s (N = 46) … {TMAX} je Übergang
(p18 brauchte bei N = 43/44 schon 122/60 Segmente). Mit p18 zusammen: **Zweig p₁ → p_{NN} rigoros als σ-Homotopie-Zweig** (Box-Radien ·ρ^k ≤ 1e-120;
Daten p24_zustand.json, Protokolle p24_lauf*.log). Weiter ist nur Rechenzeit (~M³).
Vorstudien: p22 (eine Kugel in sup über D̄ ∪ P(D̄) trägt NICHT: max_k sup_{{P(D̄)}}|J⁻¹x^k|/ρ^k = 88 → 1.6e4 für N = 20 → 40, ~1.3^N), p23 (die Richtungs-
Ableitung A = J⁻¹∂_tJ(P₀ + tΔ₁) ist winzig: ‖A‖_w = 6e-6 → 4e-9 für N = 20 → 60 — das Wachstum in p18 war reines Intervall-Wrapping).
**Fazit:** (A) ist kein unabhängiger Test, sondern: [(H) mit expliziten Konstanten ab N₀] + [Zertifikat bis N₀ — **erledigt für jedes N₀ ≤ {NN}**] +
[zertifizierte f_K-Werte auf D̄ auf ≈ R_{{N₀}} (N₀ = 100: ≈ 6e-8·C; N₀ = 150: ≈ 8e-9·C) — **fehlt, Engpass; Weg §9.8**].""")

abschnitt_98 = r"""

### 9.8 Zertifizierte Kneser-Werte (Auftrag Runde 3.1; P) — Befund und Weg

**Befund gp-tetration (nur gelesen, nichts geändert).** Dort gibt es **keine** zertifizierten Kneser-Werte zur Basis e:
- Paper VI rekonstruiert T_d aus einem als bekannt *vorausgesetzten* Anker T_e. Das „Anker-Zertifikat“ (Def. 11.3: Intervall-Tschebyschow-Darstellung
  von T_e auf einer Umgebung von [0,1] mit Restschranke) ist Voraussetzung, nicht Ergebnis; vorhanden ist nur ein Vorzertifikat mit gesetzten Puffern
  (App. A.1: Restschranke 3.8e-44, „precertificates rather than final interval certificates“). Das obere End-Kleben bei e ist ein Experiment
  (App. A.3: negative Hardy-Moden 6.9e-7, Überlapp-Residuum 5.8e-3; „evidence …, not yet the interval certificate“).
- Das 8r1-Zertifikat (PASS: Y = 3.15e-9, q = 0.150, p(r) < 0) betrifft die Basis-Strömung für b ≈ 2.043…2.060 um die Engine-Trajektorie, nicht b = e.
- research/tools/m54_cert_proto.py setzt Engine-Werte als Bälle mit ENGINE_EPS = 1e-33 an (Behauptung, kein Beweis).
- **Brauchbar:** (i) die Identifikation **Paulsen–Cowgill Prop. 2** in der genauen Fassung von Paper VI Thm 11.11: F₁, F₂ holomorph auf H_{−2} = {Re z > −2},
  F(z+1) = e^{F(z)}, F(0) = 1, für jedes x > −2 F(x+iy) → L (y → +∞) und → L̄ (y → −∞) ⇒ F₁ = F₂ (beide Grenzwerte sind Voraussetzung!);
  (ii) der Bauplan Paper VI Def. 11.8/11.12, Thm 11.13 (End-Kleben über Hardy-θ mit sup|θ'| < 1, rekursive Schnittebenen-Fortsetzung);
  (iii) Bausteine certification8r1/tmodel.py (rigorose Taylor-Modelle, Cauchy-Schwänze für exp/log, Zweigtreue), koenigs.py (Koenigs-Koordinate mit
  hergeleitetem Stopp-Rest), tails.py.

**Warum es keinen billigen Weg gibt.** (1) Ein Residuum (q∘q ≈ eˣ) bestimmt f_K nicht — jede Abel-periodische Störung ist auch Lösung; die Auswahl sitzt
allein in Realität + Grenzwerten L, L̄. (2) Knesers θ (sexp = P(L^{z+θ(z)}) oben, P = Poincaré-Funktion an L) ist an den ganzen Zahlen der reellen
Achse singulär (ln-ln-Kaskade, Paper V Lemma 5): Fourier-θ konvergiert auf ℝ nur algebraisch ⇒ man braucht zwei Karten (Taylor nahe ℝ, θ oben).
(3) Kouznetsov-Fixpunkte (Cauchy-Integral über ein Rechteck) erzeugen nicht automatisch exakte Lösungen (Randwerte ≠ Daten) ⇒ Identifikation unklar.

**Satz P-V (statisches Kneser-Zertifikat, a posteriori; Bauplan mit Beweisskizze, NICHT ausgeführt).**
*Daten:* T̃ Polynom (Taylor von sexp um 0, reelle Koeffizienten, T̃(0) = 1); θ̃(z) = Σ_{k≤K} c̃_k e^{2πikz}; Ũ(z) := P(L^{z+θ̃(z)}). **Ũ ist ganz und erfüllt
Ũ(z+1) = e^{Ũ(z)} exakt** (P(Lw) = e^{P(w)}, θ̃ periodisch). *Defekte* (Intervallrechnung, skalare Kantorovich-Schranken): d auf der Naht
{Re z = −1/2, |Im z| ≤ H+δ} mit T̃(ζ + 1 + d(ζ)) = e^{T̃(ζ)}; e auf dem Band {|Im z − H| ≤ δ} mit Ũ(ζ + e(ζ)) = T̃(ζ); η := max(sup|d|, sup|e|).
*Ansatz:* F = Ũ∘(id + v) oben (v 1-periodisch), F = T̃∘(id + w) im Rechteck |Re z| ≤ 1/2 + δ, |Im z| ≤ H + δ, unten gespiegelt. Bedingungen:
v − w = e(·+w) auf dem Band, w(z+1) − w(z) = d(z + w(z)) auf der Naht — ein **additives Cousin-Problem auf dem Zylinder ℂ/ℤ** mit Sprüngen kompakten
Trägers (Naht nur für |Im z| ≤ H, da Ũ die Funktionalgleichung exakt erfüllt; Schleifen Im z = ±H). Lösung: Ω = C + (1/2i)∫_Γ J(t)cot π(t − z) dt,
beschränkt bei ±i∞ (Liouville auf ℂ*: homogene Lösungen sind Konstanten, festgelegt durch F(0) = 1). Auswertung auf um δ/2 verschobenen Konturen ⇒
sup-Normen mit expliziter Kernkonstante C_Γ ~ (Längen)·(1/(πδ) + 1); Nichtlinearität mit Lipschitz-Konstante ~ C_Γ·sup|d'|, |e'| ≪ 1 ⇒ Kontraktion,
|v|, |w| ≤ 2C_Γη. Die Ecken-Verträglichkeit am Tripelpunkt (Naht ∩ Schleife) gilt für die nichtlineare Formulierung identisch (Kozyklus-Identität, da
Ũ exakt ist) — aber nur für eine Iterierte, die die Naht-Relation selbst exakt erfüllt; im Iterationsschritt ist die Kozyklus-Bedingung auf dem
Eck-Quadrat (Band ∩ Nahtstreifen) um O(Lip·Residuum) verletzt, und der Zylinder-Cauchy-Kern erzeugt dort einen log-Term. Das ist der eine echt
technische Punkt (Ausweg: in Band ∩ Nahtstreifen die Naht-Relation NICHT separat stellen — sie folgt dort aus Band-Relation + Periodizität von v —,
oder Newton–Kantorovich mit Projektion auf konsistente Daten plus nichtlinearer Kozyklus-Identität).
*Identifikation:* F ist holomorph auf dem Streifen, F(z+1) = e^{F(z)}, reell auf ℝ, F → L/L̄ bei ±i∞; Fortsetzung nach links über den Hauptzweig des log
verlangt F ∉ (−∞, 0] auf −1 < Re z ≤ 1/2 (bzw. F ∉ (0, 1] auf 0 < Re z ≤ 1/2 abseits ℝ) — endliche Intervallprüfung wie Paper VI Def. 11.12 ⇒ F holomorph
auf H_{−2} mit beiden Grenzwerten ⇒ **F = sexp_K** (Paulsen–Cowgill Prop. 2). *Folgerung:* |sexp_K − T̃| ≤ sup|T̃'|·2C_Γη im Rechteck, |sexp_K − Ũ| ≤
sup|Ũ'|·2C_Γη oben; f_K = sexp(slog + 1/2) auf D̄ per Intervall-Newton für slog, nahe L über Ψ = P⁻¹ (koenigs-artige Rekursion mit Stopp-Rest) und θ —
Fehler ≈ (1 + max|f_K'|)·(Kondition von slog)·2C_Γη, also wenige Vielfache von η.
*Erreichbar:* Mit Motordaten (40–50 Stellen) η ≈ 1e-30 und darunter — weit unter dem Bedarf (≈ R_{N₀}, 1e-8…1e-9). Für Os Wunsch (Koeffizienten
bis k ≈ 600 mit 120 Stellen) bräuchte man η ≈ 1e-135 (Cauchy-Formel auf |x| = 1.3 verliert (ρ/1.3)^k ≈ 1.06^k) — dafür Motordaten mit ≈ 140 Stellen.
*Baustein K1 erledigt (p25, rigoros):* zertifizierte Poincaré-Funktion: p_n (n ≤ 80) per Rekursion in arb, Rest R = P − T₈₀ als Fixpunkt von
𝒯[R] = L·e^{T(w/L)−L+R(w/L)} − T (Kontraktion, Lipschitz 1.3e-11) ⇒ |P − T₈₀| ≤ 3.3e-51·|w|^{81} auf |w| ≤ 0.5 (≤ 1.3e-75), überall per P(w) = exp^k(P(w/L^k));
Funktionalgleichung an Stichproben (|w| bis 3) auf ≤ 1e-75 eingeschlossen.
*Aufwand (ehrlich), Rest:* Zylinder-Cousin mit Konstanten, Kontraktion, Eckterm (1–2 Tage) · Schnittebene (½ Tag) ·
f_K auf D̄ einschliesslich der Umgebung von L (1 Tag): **≈ 3 Tage, eigener Auftrag**. Bis dahin sind (A), o22 und Os P1–P4 nur bedingt (Kneser-Daten
als Bälle mit geschätztem Fehler).
"""
marke = "\n---\n"
# 9.8 an das Ende von §9 (= Dateiende oder vor der nächsten Hauptüberschrift nach §9.7)
idx = s.index("### 9.7 Startbedingung (A)")
nxt = s.find("\n## ", idx)
if nxt == -1:
    s = s.rstrip("\n") + abschnitt_98
else:
    s = s[:nxt] + abschnitt_98 + s[nxt:]
p.write_text(s, encoding="utf-8")
print("ARCHITEKTUR §9 aktualisiert, N =", NN)

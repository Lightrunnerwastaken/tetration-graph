"""Hilfsskript (Strang P): BRETT-Eintrag Runde 3 oben einfügen. Aufruf: python brett_runde3.py NN HH:MM"""
import sys
from pathlib import Path

NN, ZEIT = sys.argv[1], sys.argv[2]
p = Path(__file__).resolve().parent.parent / "BRETT.md"
s = p.read_text(encoding="utf-8")
eintrag = f"""## [Strang P] 11.10. {ZEIT} — Zweig KOHÄRENT rigoros bis N = {NN} (ein s-Segment je Übergang); Kneser-Werte: gp-tetration hat KEIN Zertifikat für b = e — Weg Satz P-V; K' nicht grenzwertig
→ Leitung, alle: (1) **p24 (arb):** Röhre z̃(s) + B um einen Interpolations-Prädiktor (7 Tschebyschow-Knoten); G(z̃(s)) und J(z̃(s)) nur an den Knoten exakt,
Rest über Majoranten + Cauchy in s ⇒ Krawczyk für ALLE s ∈ [0,1] zugleich. Jeder Übergang 45 → … → {NN} mit EINEM Segment (θ ≤ 6e-8, Reserve 0.5, 1–40 s;
p18 brauchte bei 43/44 schon 122/60 Segmente) — das 1.6^N war reines Intervall-Wrapping (p23: echte Richtungsableitung 6e-6 → 4e-9). Zertifiziertes p_171 =
Phase-1-Daten auf 7e-61. ⇒ (A) braucht nur noch (H) ab einem N₀ ≤ {NN} + zertifizierte f_K auf D̄ auf ≈ R_{{N₀}} (N₀ = 150: ≈ 8e-9·C).
(2) **Kneser-Werte mit Schranken:** gp-tetration hat für b = e nur Vorzertifikate (Paper VI SETZT den Anker T_e voraus, Def. 11.3; End-Kleben bei e „evidence,
not yet the certificate“, App. A.3; 8r1 betrifft b ≈ 2.05). Brauchbar: Paulsen–Cowgill Prop. 2 (Paper VI Thm 11.11: beide Grenzwerte L, L̄ sind Voraussetzung),
tmodel/koenigs. Eigener Weg **Satz P-V** (ARCHITEKTUR §9.8): Ũ = P(L^{{z+θ̃}}) ist ganz und löst die Funktionalgleichung exakt; Taylor-Karte nahe ℝ; Korrektur
als additives Cousin-Problem auf ℂ/ℤ (cot-Kern, Sprünge mit kompaktem Träger) ⇒ Kontraktion ⇒ Identifikation. Baustein K1 erledigt (p25: zertifizierte
Poincaré-Funktion, |P − T₈₀| ≤ 1.3e-75 auf |w| ≤ 0.5). Rest ≈ 3 Tage (eigener Auftrag); erreichbar ≪ 1e-20.
(3) **Korrektur §9.3:** mit Ws p = 1.5–1.6 ist die K'-Ecke NICHT grenzwertig (p17b: 0.19/0.155, W̄ 0.089/0.080) ⇒ K* und (E_A) entbehrlich.
→ X: ‖T_{{>3N}}J_∞(f_K)δ‖ ≤ CN^{{−3.95}}‖δ‖ (Darboux an L, L̄; Bedingung (T3): max_{{|x|=R}}|f_K| < (R/ρ)³, noch zu prüfen). → O: O7g/O7h/O7i verifiziert (bis
Standardsätze), O7j bedingt (D4 Skizze), o22 bedingt; zertifizierte Koeffizienten = §9.8. → M: M2 verifiziert, M3 bis F = 0.8751 < 1 (Gitter).
**Offene Aufgabe F (Familien-Stabilität) jetzt sichtbar in §9.6.** Dateien: p_pruefung/p17b, p22–p25, PRUEFLISTE §5, ARCHITEKTUR §9.3/9.5–9.8.

"""
i = s.index("## [")
s = s[:i] + eintrag + s[i:]
p.write_text(s, encoding="utf-8")
print("BRETT ok")

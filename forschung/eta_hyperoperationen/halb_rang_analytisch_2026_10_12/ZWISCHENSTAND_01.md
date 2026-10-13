# Zwischenstand: analytische Klasse und intrinsischer quadratischer Rangjet

Stand 2026-10-12. Der H₀-Strang bleibt pausiert; dessen Dateien wurden nicht geändert.

**Bewiesen (Herleitung wird in MATHEMATIK.md ausgearbeitet):** Eine Skala holomorpher Funktionen auf rechten Halbebenen mit polynomialem Wachstum in vertikaler Richtung erlaubt eine eindeutig normierte Abel–Plana-Summation. Differentiation und Produkte bleiben mit explizitem Verlust von Halbebenenbreite bzw. Wachstumsgewichten in der Skala. Die lineare Ranggruppe auf Gamma-Moden besitzt eine analytische Integralform.

**Neue exakte Rechnung:** Für v(x)=x(x−1) sind die quadratischen Werte am ganzzahligen Höhenknoten k ohne Polynomableitung zugänglich: q_s(k)=U_k(s)(ψ(s+3)+γ)+V_k(s). U_k,V_k wurden für k=3,…,128 als rationale Polynome berechnet und ihre Differenzgleichungen koeffizientenweise geprüft. Der längste Lauf brauchte 265.50 Sekunden.

**Numerisch:** Die anschließende Newton-Fortsetzung in der Höhe gibt bei x=0.5 mit 32/64/128 Koeffizienten 0.03676976 / 0.03689300 / 0.03692866. Die alte endliche Polynomrekonstruktion hatte dort einen negativen Drift. Eine unabhängig berechnete Newton-Fortsetzung zuerst im Rang nähert sich ebenfalls diesen Werten. Das ist kein allgemeiner Satz über das Vertauschen zweier Grenzprozesse.

**Neue geschlossene Struktur:** Die erzeugende Funktion der U-Koeffizienten ist

    2 w² (1+w)^s [P_(s+2)(1+2w)−1].

P ist die Legendre-Funktion erster Art. Für die vollständige quadratische Korrektur gibt es eine Darstellung durch Jacobi-Parameterableitungen und eine normierte Rangsumme. Die Integralform wurde bei w=0.1 und 0.3 gegen 128 exakte Koeffizienten geprüft, Fehler <3·10^−41 bei 40 Dezimalstellen.

**Noch offen:** ein sauberer gleichmäßiger Randbeweis für die volle quadratische Höhenfortsetzung; die analytische nichtlineare Ranggruppe; deren Anwendung auf b^x. Die globale Klasse endlichen horizontalen Exponentialtyps kann bei b>exp(1/e) die Tetration nicht enthalten. Sie ist deshalb eine Arbeitsklasse für Rangjets und kein bereits fertiger Definitionsraum aller Hyperoperationen.

# Zwischenstand 01 — Etappe 1 (13.10.2026)

- Exakter Rechenkern `exppoly.py`: Rangschritt T als eps-Jet um x+1, alle Ordnungen, im Raum der Exponentialpolynome Σ 2^{mx}P_m(x) mit Ballarithmetik (python-flint). Geprüft gegen die direkte Rekursion an ganzen Höhen (Rang 1, 6 Ordnungen, 12 Stellen) und über die Superfunktionsgleichung an den Höhen 1…6 für alle Ränge n≤24, Ordnungen ≤10 (relative Abweichung 5·10⁻⁶⁰ = Speichergenauigkeit).
- Rang 1, eps=1, 52 Ordnungen: die Reihe **divergiert**. Glieder fallen bis ~3·10⁻⁵ (Ordnung 6–10) und wachsen dann faktoriell (log₁₀-Zuwachs pro Ordnung 0.43 bei j≈25, 0.58 bei j≈35). Die Partialsummen 4…12 schliessen Kneser an allen reellen Testpunkten ein.
- Halbrang über Newtonreihe im Rang + Levin-u: Ordnung 1 trifft die geschlossene Formel (A18) auf 34 Stellen (N=40). Ordnungen 2–7 stabil auf ~10⁻⁷…10⁻⁹.
- **Unerwartet:** Ab Ordnung 8 wachsen die ganzzahligen Rangjets c_n^{(j)}(x) exponentiell in n, mit von x unabhängigen Raten (j=10: ≈117, j=12: ≈−2667). Die Rangfortsetzung dieser Ordnungen ist mit der Newton-/Levin-Auswahl nicht bestimmt.

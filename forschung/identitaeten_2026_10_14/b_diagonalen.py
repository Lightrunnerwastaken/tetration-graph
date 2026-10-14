"""Prüft die Diagonalen der Binomialbasis-Dreiecke (Iteration von eˣ−1 und x·eˣ) für n ≤ 16:
  eˣ−1 : B(n, n−1) = n!(n−1)!/2^{n−1} (maximale Ketten im Partitionsverband), n!c_n(−1) = (−1)^{n−1}(n−1)! (Möbius)
  x·eˣ : B(n, n−1) = n!(n−1)!
"""
import json
from math import comb, factorial
from pathlib import Path
H = Path(__file__).resolve().parent
d = json.loads((H / "daten" / "b_kombinatorik.json").read_text(encoding="utf-8"))["daten"]
e, x = d["expm1"]["dreieck_binomialbasis"], d["xexp"]["dreieck_binomialbasis"]
N = len(e)
print("expm1 Diagonale:", all(e[n - 1][n - 1] == factorial(n) * factorial(n - 1) // 2 ** (n - 1) for n in range(2, N + 1)))
print("expm1 t=-1 (Möbius):", all(sum(b * (-1) ** j for j, b in enumerate(e[n - 1])) == (-1) ** (n - 1) * factorial(n - 1) for n in range(1, N + 1)))
print("xexp Diagonale:", all(x[n - 1][n - 1] == factorial(n) * factorial(n - 1) for n in range(2, N + 1)))
# x·eˣ bei t = −1: Umkehrfunktion W(x): n![x^n]W = (−n)^{n−1}
print("xexp t=-1 (Lambert W):", all(sum(b * (-1) ** j for j, b in enumerate(x[n - 1])) == (-n) ** (n - 1) for n in range(1, N + 1)))

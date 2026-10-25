"""Gebiet 02, Teil F: Modul. |n a_n| gegen |h'(Re c − ln n + i0)| = |1+θ'(σ_n)|/sexp'(σ_n−2)  (Basis e oder 2).
Aufruf: python p4_modul.py [e|2]"""
import sys, json
import mpmath as mp
import p2_formeln as P2   # lädt Daten, h', c (gibt dabei die p2-Tabellen aus)
print("\n== Modul ==")
print("    n | n|a_n|      | |h'(v)| führend | |h'(Re v+i0)| Achse | Achse/führend")
for n in [16, 256, 1024, 4096]:
    v = (P2.c - mp.log(n)) / P2.lb
    f0 = abs(P2.hp(v) / P2.lb)
    f1 = abs(P2.hp(mp.mpc(mp.re(v), P2.EPS)) / P2.lb)
    print(f"{n:5d} | {float(n*abs(P2.a[n])):.6e} | {float(f0):.6e} | {float(f1):.6e} | {float(f1/f0):.5f}")

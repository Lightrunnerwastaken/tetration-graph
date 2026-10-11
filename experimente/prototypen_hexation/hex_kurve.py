"""Kurve γ(w) = α_q(w) (F_q(γ) = w) auf einem gleichmässigen w-Raster; |γ'(w)| = 1/|F_q'|. Startwerte aus der groben Kurve (K = 128)."""
import sys, time, json
sys.path.insert(0, r"C:/Users/janis/AppData/Local/Temp/claude/D--Users-janis-Documents-Matura/16e6767d-2c78-46bf-8ec0-6069646a66c1/scratchpad")
import mpmath as mp
from hexq import HexQ
t0 = time.time()
def log(*a):
    print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)

H = HexQ("2", mp.mpf("-1.7439"), mp.mpc("-1.99294843593232", "1.44512849663921"))
alt = json.load(open(r"C:/Users/janis/AppData/Local/Temp/claude/D--Users-janis-Documents-Matura/16e6767d-2c78-46bf-8ec0-6069646a66c1/scratchpad/hex_lauf_2.json"))
galt = [mp.mpf(x) for x in alt["128"]["g"]]
talt = [mp.mpf(j) / 128 for j in range(128)]
Falt = [mp.re(v) for v in H.F([t + 1j * g for t, g in zip(talt, galt)])]
log("F entlang der alten Kurve monoton:", all(Falt[i] < Falt[i + 1] for i in range(127)))
# Wo ist die alte Kurve steil? Abstand aufeinanderfolgender Punkte gegen Wertzuwachs
for j in range(0, 128, 4):
    dz = abs((talt[j + 1] - talt[j]) + 1j * (galt[(j + 1) % 128] - galt[j])) if j < 127 else 0
    dw = Falt[j + 1] - Falt[j] if j < 127 else 0
    log(f"t = {mp.nstr(talt[j], 4)}  g = {mp.nstr(galt[j], 8)}  F = {mp.nstr(Falt[j], 8)}  |dζ|/dw = {mp.nstr(dz / dw if dw else 0, 4)}")
h = mp.mpf("1e-12")
ableit = H.F([t + 1j * g + h for t, g in zip(talt, galt)])
dF = [(a - b) / h for a, b in zip(ableit, H.F([t + 1j * g for t, g in zip(talt, galt)]))]
jmin = min(range(128), key=lambda j: abs(dF[j]))
log("kleinstes |F'| auf der alten Kurve:", mp.nstr(abs(dF[jmin]), 4), "bei t =", mp.nstr(talt[jmin], 4), "F =", mp.nstr(Falt[jmin], 8))
log("|F'| entlang:", [mp.nstr(abs(d), 3) for d in dF[::8]])

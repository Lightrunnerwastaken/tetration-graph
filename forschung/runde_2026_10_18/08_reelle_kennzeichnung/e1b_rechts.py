"""E1b: Kneser rechts von 0.5: exp-Weg über zwei verschiedene Zentren; verlässliche Ordnung = wo beide auf 1e-6 übereinstimmen."""
import mpmath as mp, rk_lib as R
R.init(60); N = 100
def k_von(x, zweit):
    y = x - 1
    cs = sorted(R._Z, key=lambda c: abs(y - c))
    c0 = cs[1] if zweit else cs[0]
    return R.exp_s(R.shift(R._Z[c0], y - c0, N), N)
x = mp.mpf('0.55')
while x <= 1.5001:
    a = k_von(x, False); b = k_von(x, True)
    K = next((k for k in range(1, N) if abs(a[k] - b[k]) > mp.mpf('1e-6') * abs(a[k])), N)
    neg = [k for k in range(1, K, 2) if a[k] <= 0]
    print(mp.nstr(x, 4), "verlässlich bis k <", K, " ungerade <= 0:", neg, " gerade>0 bis", max([k for k in range(2, K, 2) if a[k] > 0] or [0]), flush=True)
    x += mp.mpf('0.05')

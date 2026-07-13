# wildlamp/quant/dye_pr.py

import math

PR_BASELINE = (220, 80, 120)

def color_pr_from_copies(N: float):
    pH0 = 8.3
    alpha = 0.12
    pH = pH0 - alpha * math.log10(max(N, 1))
    x = 1 / (1 + math.exp(-(pH - 7.4) * 4))
    r0, g0, b0 = PR_BASELINE
    r1, g1, b1 = 240, 220, 60
    return (
        int(r0 + (1 - x) * (r1 - r0)),
        int(g0 + (1 - x) * (g1 - g0)),
        int(b0 + (1 - x) * (b1 - b0))
    )

def delta_rgb(R, G, B, R0, G0, B0):
    return abs(R - R0) + abs(G - G0) + abs(B - B0)

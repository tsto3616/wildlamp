import numpy as np
from wildlamp.simulations.model_rtlamp import simulate_rt_qlamp
from wildlamp.simulations.dye import color_pr_from_copies, delta_rgb, PR_BASELINE

QUANT_CONC = [1e6, 3e5, 1e5, 3e4, 1e4, 3e3, 1e3, 3e2, 1e2]

def simulate_quant_curve_pr(FIP, BIP, outer, template_seq):
    logN_list = []
    tt_list = []
    tt_colour_list = []

    R0, G0, B0 = PR_BASELINE

    for N0 in QUANT_CONC:
        times, copies, tt = simulate_rt_qlamp(FIP, BIP, outer, template_seq, N0=N0)

        tt_colour = None
        for t, N in zip(times, copies):
            R, G, B = color_pr_from_copies(N)
            if delta_rgb(R, G, B, R0, G0, B0) >= 40:
                tt_colour = t
                break

        logN_list.append(np.log10(N0))
        tt_list.append(tt)
        tt_colour_list.append(tt_colour)

    return np.array(logN_list), np.array(tt_list), np.array(tt_colour_list)

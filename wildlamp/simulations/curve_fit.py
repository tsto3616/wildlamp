# wildlamp/quant/curve_fit.py

import numpy as np

def fit_standard_curve(logN, tt):
    tt_arr = np.array([np.nan if v is None else v for v in tt], dtype=float)
    mask = ~np.isnan(tt_arr)
    slope, intercept = np.polyfit(logN[mask], tt_arr[mask], 1)
    return intercept, slope

def copies_from_Tt(tt, intercept, slope):
    logN = (tt - intercept) / slope
    return 10 ** logN

"""
Thermodynamic scoring utilities for primer optimisation.

Provides:
    thermo_score(seq) → (tm, score)
"""

import primer3
from wildlamp.utils.seq import gc  


def thermo_score(seq):
    """
    Compute thermodynamic score for a primer sequence.

    Returns
    -------
    tm : float
        Melting temperature from primer3.
    score : float
        Composite thermodynamic score:
            tm
          - hairpin_tm
          - dimer_tm
          - GC penalty (outside 35–70%)
    """

    ta = primer3.thermoanalysis.ThermoAnalysis()

    tm = ta.calc_tm(seq)
    hp = ta.calc_hairpin(seq)
    hd = ta.calc_homodimer(seq)

    hairpin_tm = hp.tm if hp.structure_found else 0.0
    dimer_tm   = hd.tm if hd.structure_found else 0.0

    g = gc(seq)

    # GC penalty outside 35–70%
    gc_pen = 0.0
    if g < 0.35:
        gc_pen = (0.35 - g) * 100
    elif g > 0.70:
        gc_pen = (g - 0.70) * 100

    score = tm - hairpin_tm - dimer_tm - gc_pen

    return tm, score

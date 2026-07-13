"""
Primer3 thermodynamic evaluation utilities for the wildlamp package.

Provides:
- evaluate_primer3(): compute TM, hairpin TM, and homodimer TM for each
  primer in a LAMP set using primer3-py bindings.
"""

import primer3


def evaluate_primer3(s: dict) -> dict:
    """
    Evaluate thermodynamic properties of primers in a LAMP set using Primer3.

    Parameters
    ----------
    s : dict
        Dictionary containing primer sequences under keys:
        "F3", "F2", "F1", "B3", "B2", "B1".

    Returns
    -------
    dict
        {
            primer_name: {
                "tm": float,
                "hairpin_tm": float,
                "dimer_tm": float
            }
        }
    """

    ta = primer3.thermoanalysis.ThermoAnalysis()
    out = {}

    for r in ["F3", "F2", "F1", "B3", "B2", "B1"]:
        p = s[r]

        tm = ta.calc_tm(p)
        hp = ta.calc_hairpin(p)
        hd = ta.calc_homodimer(p)

        out[r] = {
            "tm": tm,
            "hairpin_tm": hp.tm if hp.structure_found else 0.0,
            "dimer_tm": hd.tm if hd.structure_found else 0.0
        }

    return out

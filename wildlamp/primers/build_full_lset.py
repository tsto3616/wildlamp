"""
Build a complete LAMP primer set (lset) from original + optimised primers.

Used by final_primers() to assemble:
- F3, F2, F1, B3, B2, B1
- Their genomic positions
- FIP and BIP constructed from optimised primers
"""

from wildlamp.utils.seq import rc


def build_full_lset(primers, optimised, ref):
    """
    Construct a full LAMP primer set dictionary.

    Parameters
    ----------
    primers : dict
        Original primer sequences:
            {"F3":..., "F2":..., "F1":..., "B3":..., "B2":..., "B1":...}

    optimised : dict
        Optimised primer info:
            {
                "F3": {"seq":..., "pos":..., ...},
                ...
            }

    ref : str
        Reference sequence (cleaned alignment first sequence)

    Returns
    -------
    dict
        {
            "F3": seq, "F3_pos": pos,
            "F2": seq, "F2_pos": pos,
            "F1": seq, "F1_pos": pos,
            "B3": seq, "B3_pos": pos,
            "B2": seq, "B2_pos": pos,
            "B1": seq, "B1_pos": pos,
            "FIP": rc(F1) + F2,
            "BIP": rc(B1) + B2
        }
    """

    lset = {}

    for r in ["F3", "F2", "F1", "B3", "B2", "B1"]:
        # Use optimised primer if available, otherwise original
        seq = optimised[r]["seq"] if r in optimised else primers[r]

        # Use optimised position if available, otherwise find original in ref
        pos = optimised[r]["pos"] if r in optimised else ref.find(primers[r])

        lset[r] = seq
        lset[f"{r}_pos"] = pos

    # Build FIP and BIP from optimised sequences
    lset["FIP"] = rc(lset["F1"]) + lset["F2"]
    lset["BIP"] = rc(lset["B1"]) + lset["B2"]

    return lset

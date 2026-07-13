"""
Amplicon heuristic and amplicon scoring for the wildlamp package.

Provides:
- amplicon_heuristic(): lightweight amplicon quality scoring based on
  length, GC content, and homopolymer penalties.
- amplicon_score(): convert amplicon_heuristic() output into a numeric score.
"""

from utils.seq import gc, bad_run


def amplicon_heuristic(seq: str):
    """
    Compute a simple amplicon quality heuristic.

    Parameters
    ----------
    seq : str
        Amplicon DNA sequence.

    Returns
    -------
    dict
        {
            "length": int,
            "gc": float,
            "penalty": float
        }

    Notes
    -----
    Penalty components:
    - length_pen: deviation from ideal ~200 bp
    - gc_pen: penalty for GC outside 40–65%
    - homo_pen: penalty for homopolymers >6 bases
    """

    L = len(seq)
    g = gc(seq)

    # Length penalty: ideal ~200 bp, acceptable 120–280 bp
    length_pen = abs(L - 200) * 0.5 if not (120 <= L <= 280) else 0

    # GC penalty
    gc_pen = 0
    if g < 0.40:
        gc_pen = (0.40 - g) * 100
    if g > 0.65:
        gc_pen = (g - 0.65) * 100

    # Homopolymer penalty
    homo_pen = 30 if bad_run(seq) else 0

    return {
        "length": L,
        "gc": g,
        "penalty": length_pen + gc_pen + homo_pen
    }

def amplicon_score(amp_info: dict) -> float:
    """
    Compute a numeric amplicon score from heuristic components.

    Parameters
    ----------
    amp_info : dict
        Output of amplicon_heuristic(), containing:
        - "length": int
        - "gc": float
        - "penalty": float

    Returns
    -------
    float
        Amplicon score combining length suitability, GC suitability,
        and penalty deductions.
    """

    L = amp_info["length"]
    g = amp_info["gc"]
    pen = amp_info["penalty"]

    score = 0

    # Length suitability
    if 140 <= L <= 240:
        score += 40
    elif 100 <= L <= 280:
        score += 20

    # GC suitability
    if 0.40 <= g <= 0.65:
        score += 20

    # Penalty deduction
    score -= 0.5 * pen

    return score

"""
Mismatch scoring utilities for the wildlamp package.

Provides:
- mismatch_cost(): weighted mismatch penalty between primer and target segment.
"""

def mismatch_cost(primer: str, seg: str) -> float:
    """
    Compute mismatch penalty between a primer and a target segment.

    Parameters
    ----------
    primer : str
        Primer sequence.
    seg : str
        Target sequence segment of equal length.

    Returns
    -------
    float
        Weighted mismatch penalty.
    """
    cost = 0.0
    L = len(primer)

    for i, (a, b) in enumerate(zip(primer, seg), start=1):

        if a == b:
            continue

        if b == "-":
            c_type = 3
        elif (a, b) in [("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")]:
            c_type = 1
        else:
            c_type = 2

        w_pos = 2.0 if i > L - 5 else 1.0
        cost += w_pos * c_type

    return cost

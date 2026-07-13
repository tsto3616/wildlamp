"""
Degenerate-base suggestion utilities for the wildlamp package.

Provides:
- DEGEN_MAP: mapping from sets of bases → IUPAC degenerate code
- suggest_degenerate(): replace a primer position with appropriate IUPAC code
"""

DEGEN_MAP = {
    frozenset(["A","G"]): "R",
    frozenset(["C","T"]): "Y",
    frozenset(["G","C"]): "S",
    frozenset(["A","T"]): "W",
    frozenset(["G","T"]): "K",
    frozenset(["A","C"]): "M",
    frozenset(["A","C","G"]): "V",
    frozenset(["A","C","T"]): "H",
    frozenset(["A","G","T"]): "D",
    frozenset(["C","G","T"]): "B",
    frozenset(["A","C","G","T"]): "N",
}


def suggest_degenerate(primer: str, idx: int, bases: list[str]) -> str:
    """
    Suggest a degenerate-base replacement at a mismatch position.

    Parameters
    ----------
    primer : str
        Original primer sequence.
    idx : int
        Index of mismatch within primer.
    bases : list[str]
        Bases observed across species at this genomic position.

    Returns
    -------
    str
        Primer with IUPAC degenerate base inserted, or original primer.
    """
    base_set = frozenset(bases)
    if base_set in DEGEN_MAP:
        deg = DEGEN_MAP[base_set]
        return primer[:idx] + deg + primer[idx+1:]
    return primer

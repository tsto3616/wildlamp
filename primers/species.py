"""
Species‑coverage utilities for the wildlamp package.

Provides:
- primer_species_coverage(): species that a primer can amplify
- filter_primers_by_species(): retain primers amplifying ≥ min_species
"""

from primers.mismatch import mismatch_cost


def primer_species_coverage(
    primer: str,
    pos: int,
    aln: list[str],
    names: list[str],
    cost_max: float = 4.0
):
    """
    Determine which species a primer can amplify based on mismatch cost.

    Parameters
    ----------
    primer : str
        Primer sequence.
    pos : int
        Primer binding position in the reference sequence.
    aln : list[str]
        Cleaned alignment sequences.
    names : list[str]
        Species names corresponding to aln.
    cost_max : float
        Maximum allowed mismatch cost.

    Returns
    -------
    list[str]
        Species names for which the primer is acceptable.
    """
    ok = []

    for i, seq in enumerate(aln):
        seg = seq[pos:pos + len(primer)]
        if len(seg) != len(primer):
            continue
        if mismatch_cost(primer, seg) <= cost_max:
            ok.append(names[i])

    return ok


def filter_primers_by_species(
    primers: list[tuple[int, str]],
    aln: list[str],
    names: list[str],
    cost_max: float = 4.0,
    min_species: int = 5
):
    """
    Filter primer candidates by requiring amplification across multiple species.

    Parameters
    ----------
    primers : list[(pos, primer)]
        Primer candidates.
    aln : list[str]
        Alignment sequences.
    names : list[str]
        Species names.
    cost_max : float
        Maximum mismatch cost allowed.
    min_species : int
        Minimum number of species required for retention.

    Returns
    -------
    list[(pos, primer)]
        Primers that amplify at least min_species species.
    """
    return [
        (pos, p)
        for pos, p in primers
        if len(primer_species_coverage(p, pos, aln, names, cost_max)) >= min_species
    ]

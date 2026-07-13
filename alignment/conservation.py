"""
wildlamp.alignment.conservation

This module provides conservation scoring utilities for multi-species
alignments. It is part of the wildlamp RT-LAMP primer design pipeline.

The main function, compute_conservation(), calculates how conserved each
column of an alignment is relative to a reference sequence. This is used
to prioritise primer binding sites that are stable across species.
"""

from typing import List


def compute_conservation(aln: List[str], ref: str) -> List[float]:
    """
    Compute per-column conservation relative to a reference sequence.

    Parameters
    ----------
    aln : list[str]
        A list of aligned DNA sequences (all padded to equal length).
        These sequences represent multiple species in a conservation
        genetics context.

    ref : str
        The reference sequence (usually aln[0]) used to evaluate
        conservation. Each position in ref is compared to the same
        position in all other sequences.

    Returns
    -------
    cons : list[float]
        A list of conservation scores between 0.0 and 1.0.
        cons[j] = proportion of species whose nucleotide at column j
        matches the reference nucleotide.

    Notes
    -----
    - Non-ACGT characters (e.g., gaps) in the reference sequence
      automatically yield a conservation score of 0.0.
    - This function is intentionally simple and deterministic so that
      conservation scoring is reproducible for MER reviewers.
    """

    N = len(aln)      # number of species
    L = len(ref)      # alignment length

    cons = []

    for j in range(L):
        r = ref[j]

        # If the reference has a gap or ambiguous base, conservation is undefined → 0.0
        if r not in "ACGT":
            cons.append(0.0)
            continue

        # Count how many species match the reference at this position
        matches = sum(1 for seq in aln if seq[j] == r)

        # Proportion of species matching the reference
        cons.append(matches / N)

    return cons

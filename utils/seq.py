"""
Sequence utility functions for the wildlamp package.

These helpers provide:
- gc(): GC fraction
- bad_run(): detection of long homopolymer runs
- rc(): reverse complement

They are used across multiple modules:
- primer scanning
- amplicon scoring
- LAMP geometry construction
- thermodynamic filtering
"""

from typing import List


def gc(seq: str) -> float:
    """
    Compute GC fraction of a DNA sequence.

    Parameters
    ----------
    seq : str
        DNA sequence (ACGT only).

    Returns
    -------
    float
        GC fraction between 0.0 and 1.0.
    """
    return (seq.count("G") + seq.count("C")) / len(seq)


def bad_run(seq: str, max_run: int = 6) -> bool:
    """
    Detect homopolymer runs longer than max_run.

    Parameters
    ----------
    seq : str
        DNA sequence.
    max_run : int
        Maximum allowed run length.

    Returns
    -------
    bool
        True if a homopolymer run exceeds max_run.
    """
    run = 1
    for a, b in zip(seq, seq[1:]):
        if a == b:
            run += 1
            if run > max_run:
                return True
        else:
            run = 1
    return False


def rc(seq: str) -> str:
    """
    Reverse complement a DNA sequence.

    Parameters
    ----------
    seq : str
        DNA sequence.

    Returns
    -------
    str
        Reverse complement of the sequence.
    """
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

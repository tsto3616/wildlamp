"""
Primer scanning utilities for the wildlamp package.

These functions identify candidate primer windows in a reference sequence
and apply basic filtering rules (length, GC content, homopolymers, ACGT-only).
"""

from wildlamp.utils.seq import gc, bad_run


def primer_ok(p: str) -> bool:
    """
    Basic primer suitability filter.

    Parameters
    ----------
    p : str
        Primer sequence.

    Returns
    -------
    bool
        True if primer passes length, GC, nucleotide, and homopolymer checks.
    """
    if not (18 <= len(p) <= 35):
        return False
    if any(b not in "ACGT" for b in p):
        return False
    if not (0.35 <= gc(p) <= 0.70):
        return False
    if bad_run(p):
        return False
    return True


def scan_primers(ref: str):
    """
    Scan a reference sequence for candidate primers using a sliding window.

    Parameters
    ----------
    ref : str
        Reference DNA sequence.

    Returns
    -------
    list[tuple[int, str]]
        List of (position, primer_sequence) tuples.
    """
    primers = []
    pos = 0

    while pos < len(ref):
        found = False

        for k in range(18, 36):
            if pos + k > len(ref):
                break

            p = ref[pos:pos + k]

            if primer_ok(p):
                primers.append((pos, p))
                pos += max(k, 15)
                found = True
                break

        if not found:
            pos += 1

    return primers

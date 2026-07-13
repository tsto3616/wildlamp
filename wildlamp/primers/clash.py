"""
Primer clash detection utilities for the wildlamp package and primer clash-map. 

Provides:
- primers_clash(): determine whether two primers clash based on spacing rules.
- build_clash_map(): build adjacency sets of clashing primers.
"""

def primers_clash(pos1: int, len1: int, pos2: int, len2: int) -> bool:
    """
    Determine whether two primers clash based on genomic positions and lengths.

    Parameters
    ----------
    pos1 : int
        Start position of primer 1.
    len1 : int
        Length of primer 1.
    pos2 : int
        Start position of primer 2.
    len2 : int
        Length of primer 2.

    Returns
    -------
    bool
        True if primers clash (overlap or too close), False otherwise.
    """
    s1, e1 = pos1, pos1 + len1
    s2, e2 = pos2, pos2 + len2

    # Non-overlapping but too close (<15 bp)
    if e1 <= s2:
        return (s2 - e1) < 15
    if e2 <= s1:
        return (s1 - e2) < 15

    # Overlapping primers always clash
    return True

def build_clash_map(primers):
    """
    Build a clash map for a list of primers.

    Parameters
    ----------
    primers : list[(pos, seq)]
        List of primer candidates.

    Returns
    -------
    dict[int, set[int]]
        Dictionary mapping primer index → set of indices that clash with it.
    """
    clash = {}

    for i, (p1_pos, p1_seq) in enumerate(primers):
        clash[i] = set()

        for j, (p2_pos, p2_seq) in enumerate(primers):
            if i == j:
                continue

            if primers_clash(p1_pos, len(p1_seq), p2_pos, len(p2_seq)):
                clash[i].add(j)

    return clash

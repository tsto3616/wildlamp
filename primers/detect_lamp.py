"""
LAMP primer-set detection for the wildlamp package.

This module identifies valid LAMP primer sets (F3, F2, F1, B1, B2, B3)
within a reference sequence using spacing rules, clash maps, and
species-coverage scoring.
"""

from ..utils.seq import rc
from ..primers.species import primer_species_coverage


def detect_lamp_sets(ref, primers, aln, names, clash):
    """
    Detect valid LAMP primer sets from candidate primers.

    Parameters
    ----------
    ref : str
        Reference DNA sequence.
    primers : list[(pos, seq)]
        Primer candidates.
    aln : list[str]
        Cleaned alignment sequences.
    names : list[str]
        Species names.
    clash : dict[int, set[int]]
        Clash map: primer index → set of clashing primer indices.

    Returns
    -------
    list[dict]
        List of LAMP sets, each containing F3/F2/F1/B1/B2/B3 and metadata.
    """

    best = []
    best_cov = 0
    L = len(ref)

    def cov(s):
        return len(primer_species_coverage(s, aln, names))

    # Sliding windows across the reference
    for win_start in range(0, L, 100):
        win_end = win_start + 500

        window = [
            (i, pos, p)
            for i, (pos, p) in enumerate(primers)
            if win_start <= pos <= win_end
        ]

        if len(window) < 6:
            continue

        # F3
        for f3_i, f3_pos, f3 in window:

            # F2
            for f2_i, f2_pos, f2 in window:
                if f2_pos <= f3_pos:
                    continue
                if not (40 <= f2_pos - f3_pos <= 300):
                    continue
                if f2_i in clash[f3_i]:
                    continue

                # F1
                for f1_i, f1_pos, f1 in window:
                    if f1_pos <= f2_pos:
                        continue
                    if not (20 <= f1_pos - f2_pos <= 200):
                        continue
                    if f1_i in clash[f2_i] or f1_i in clash[f3_i]:
                        continue

                    F1c = rc(f1)

                    # B3
                    for b3_i, b3_pos, b3 in window:
                        if b3_pos <= f1_pos:
                            continue
                        if not (100 <= b3_pos - f1_pos <= 600):
                            continue
                        if b3_i in clash[f1_i] or b3_i in clash[f2_i] or b3_i in clash[f3_i]:
                            continue

                        # B2
                        for b2_i, b2_pos, b2 in window:
                            if b2_pos <= b3_pos:
                                continue
                            if not (20 <= b2_pos - b3_pos <= 200):
                                continue
                            if b2_i in clash[b3_i]:
                                continue

                            # B1
                            for b1_i, b1_pos, b1 in window:
                                if b1_pos <= b2_pos:
                                    continue
                                if not (20 <= b1_pos - b2_pos <= 200):
                                    continue
                                if b1_i in clash[b2_i] or b1_i in clash[b3_i]:
                                    continue

                                B1c = rc(b1)
                                FIP = F1c + f2
                                BIP = B1c + b2

                                s = {
                                    "F3": f3, "F3_pos": f3_pos,
                                    "F2": f2, "F2_pos": f2_pos,
                                    "F1": f1, "F1_pos": f1_pos, "F1c": F1c,
                                    "B3": b3, "B3_pos": b3_pos,
                                    "B2": b2, "B2_pos": b2_pos,
                                    "B1": b1, "B1_pos": b1_pos, "B1c": B1c,
                                    "FIP": FIP, "BIP": BIP,
                                    "amplicon_len": (b2_pos + len(b2)) - f2_pos
                                }

                                c = cov(s)

                                if c < best_cov:
                                    continue

                                if c > best_cov:
                                    best_cov = c
                                    best = []

                                best.append(s)

                                if len(best) > 50:
                                    best = best[-50:]

    return best

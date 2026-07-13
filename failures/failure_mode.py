"""
Failure-mode analysis for LAMP primer sets.

Identifies:
- species-specific mismatches
- mismatch positions
- alternative bases across alignment
- suggested degenerate-base corrections
"""

from ..primers.mismatch import mismatch_cost
from ..failures.suggest import suggest_degenerate


def analyze_failure_modes(lset, aln, names, cost_max=4.0):
    """
    Analyze mismatch failure modes for a LAMP primer set.

    Parameters
    ----------
    lset : dict
        LAMP set containing primers and positions.
    aln : list[str]
        Alignment sequences.
    names : list[str]
        Species names.
    cost_max : float
        Maximum allowed mismatch cost.

    Returns
    -------
    list[dict]
        Failure-mode entries describing mismatch positions and suggested fixes.
    """

    failures = []
    roles = ["F3","F2","F1","B1","B2","B3"]

    for i, seq in enumerate(aln):
        species_name = names[i]

        for r in roles:
            p = lset[r]
            pos = lset[f"{r}_pos"]
            seg = seq[pos:pos+len(p)]

            if len(seg) != len(p):
                continue

            cost = mismatch_cost(p, seg)
            if cost <= cost_max:
                continue

            mismatch_positions = [
                idx for idx, (a, b) in enumerate(zip(p, seg))
                if a != b and b in "ACGT"
            ]

            for mp in mismatch_positions:
                genomic_pos = pos + mp

                bases = []
                for sseq in aln:
                    if genomic_pos < len(sseq):
                        b = sseq[genomic_pos]
                        if b in "ACGT":
                            bases.append(b)

                if not bases:
                    continue

                alt_bases = sorted(set(bases))
                suggested = suggest_degenerate(p, mp, alt_bases)

                failures.append({
                    "species": species_name,
                    "role": r,
                    "primer_pos": pos,
                    "primer_index": mp,
                    "original_primer": p,
                    "ref_base": p[mp],
                    "alt_bases": "".join(alt_bases),
                    "suggested_primer": suggested,
                    "mismatch_cost": cost
                })

    return failures

"""
Standalone optimisation command for LAMP primers.

This function depends on three external package functions:
- primer_ok
- primer_species_coverage
- thermo_score

Users call:
    from wildlamp.optimisation.optimise_primer import optimise_primer
    best_seq, best_pos, best_score, best_cov = optimise_primer(...)

"""

from ..primers.scan import primer_ok
from ..primers.species import primer_species_coverage
from ..thermodynamics.scoring import thermo_score


def optimise_primer(
    aln, names, ref,
    original_primer,
    window_shift=10,
    tm_min=55.0,
    tm_max=62.0,
    max_mutations=2,
    no_mutate_3prime=5,
    cost_max=4.0,
    min_binding_transcripts=5
):
    """
    Optimise a primer by sliding a window across the reference sequence
    and performing limited base mutations to improve thermodynamics and
    species coverage.

    Returns:
        best_seq : str or None
        best_pos : int or None
        best_score : float or None
        best_cov : int or None
    """

    L = len(original_primer)
    pos0 = ref.find(original_primer)

    if pos0 == -1:
        print(f"[WARN] Primer {original_primer} not found in cleaned ref.")
        return None, None, None, None

    best_seq = None
    best_score = float("-inf")
    best_pos = None
    best_cov = 0

    # ------------------------------------------------------------
    # WINDOW SHIFTING
    # ------------------------------------------------------------
    for shift in range(-window_shift, window_shift + 1):
        pos = pos0 + shift

        if pos < 0 or pos + L > len(ref):
            continue

        window = ref[pos:pos + L]
        if any(b not in "ACGT" for b in window):
            continue

        base = window

        # Basic primer validity
        if not primer_ok(base, min_len=L, max_len=L):
            continue

        # Species coverage
        cov_names = primer_species_coverage(
            base, pos, aln, names, cost_max=cost_max
        )
        if len(cov_names) < min_binding_transcripts:
            continue

        # Thermodynamics
        tm, thermo = thermo_score(base)
        if not (tm_min <= tm <= tm_max):
            continue

        # Local best before mutation
        best_local_seq = base
        best_local_score = thermo + 10.0 * len(cov_names)
        best_local_cov = len(cov_names)

        # ------------------------------------------------------------
        # MUTATION LOOP
        # ------------------------------------------------------------
        for i in range(L):
            if i >= L - no_mutate_3prime:
                continue

            for b in "ACGT":
                if b == base[i]:
                    continue

                cand = base[:i] + b + base[i+1:]

                # Mutation count
                muts = sum(1 for a, c in zip(base, cand) if a != c)
                if muts > max_mutations:
                    continue

                if not primer_ok(cand, min_len=L, max_len=L):
                    continue

                cov_names_cand = primer_species_coverage(
                    cand, pos, aln, names, cost_max=cost_max
                )
                cov_count = len(cov_names_cand)
                if cov_count < min_binding_transcripts:
                    continue

                tm_c, thermo_c = thermo_score(cand)
                if not (tm_min <= tm_c <= tm_max):
                    continue

                score_c = thermo_c + 10.0 * cov_count

                if score_c > best_local_score:
                    best_local_seq = cand
                    best_local_score = score_c
                    best_local_cov = cov_count

        # ------------------------------------------------------------
        # UPDATE GLOBAL BEST
        # ------------------------------------------------------------
        if best_local_score > best_score:
            best_seq = best_local_seq
            best_score = best_local_score
            best_pos = pos
            best_cov = best_local_cov

    return best_seq, best_pos, best_score, best_cov

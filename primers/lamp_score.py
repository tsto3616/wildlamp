"""
Primer scoring utilities for the wildlamp package.

Provides:
- primer_basic_score()
- primer_thermo_score()
- primer_struct_score()
- primer_conservation_score()
- lamp_score()
- rtlamp_score()
"""

from utils.seq import gc, bad_run


def primer_basic_score(p: str) -> int:
    """
    Basic primer suitability score based on GC, length, and homopolymers.
    """
    score = 0
    if 0.40 <= gc(p) <= 0.65:
        score += 20
    if 18 <= len(p) <= 30:
        score += 20
    if not bad_run(p):
        score += 10
    return score


def primer_thermo_score(ev: dict) -> int:
    """
    Score thermodynamic suitability using Primer3 evaluation results.
    """
    tm = ev["tm"]
    hairpin_tm = ev["hairpin_tm"]
    dimer_tm = ev["dimer_tm"]

    score = 0

    if 60 <= tm <= 68:
        score += 20
    elif 55 <= tm <= 72:
        score += 10

    if hairpin_tm > 40:
        score -= 10
    if dimer_tm > 40:
        score -= 10

    return score


def primer_struct_score(vres: dict) -> int:
    """
    Score ViennaRNA structural stability.
    """
    return 10 if vres["mfe"] > -5.0 else 0


def primer_conservation_score(pos: int, length: int, cons_array: list[float]) -> float:
    """
    Score conservation of primer binding site across species.
    """
    c = sum(cons_array[pos + k] for k in range(length)) / length
    return 20 * c


def lamp_score(s: dict, cons_array: list[float], N_species: int) -> float:
    """
    Combined LAMP primer set score including:
    - basic primer scores
    - thermodynamic scores
    - structural scores
    - conservation scores
    - species coverage
    - amplicon score
    """
    roles = ["F3", "F2", "F1", "B1", "B2", "B3"]
    primer_score = 0

    for r in roles:
        p = s[r]
        pos = s[f"{r}_pos"]

        basic = primer_basic_score(p)
        thermo = primer_thermo_score(s["primer3_eval"][r])
        struct = primer_struct_score(s["vienna_validation"][r])
        cons = primer_conservation_score(pos, len(p), cons_array)

        primer_score += basic + thermo + struct + cons

    cov_score = 100.0 * (s["species_matched"] / N_species)
    amp_score_val = s["amp_score"]

    return primer_score + cov_score + amp_score_val


def rtlamp_score(s: dict) -> float:
    """
    RT‑LAMP scoring variant using final_score and amplicon penalties.
    """
    score = s["final_score"]
    L = s["amplicon_len"]

    if 140 <= L <= 240:
        score += 80
    elif 100 <= L <= 280:
        score += 40

    score -= s["amp_penalty"] * 0.5

    return score

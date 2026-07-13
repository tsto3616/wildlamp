import numpy as np
import math
import primer3.bindings as p3
from ..utils.seq import rc

def best_match_mismatches(primer: str, template: str) -> int:
    Lp = len(primer)
    Lt = len(template)
    min_mm = Lp
    for i in range(Lt - Lp + 1):
        window = template[i:i+Lp]
        mm = sum(1 for a, b in zip(primer, window) if a != b)
        if mm < min_mm:
            min_mm = mm
            if min_mm == 0:
                break
    return min_mm

def compute_hairpin(seq: str) -> float:
    return p3.calc_hairpin(seq).dg

def compute_homodimer(seq: str) -> float:
    return p3.calc_homodimer(seq).dg

def compute_heterodimer(seqA: str, seqB: str) -> float:
    return p3.calc_heterodimer(seqA, seqB).dg

def has_3prime_dimer(seqA: str, seqB: str, min_match: int = 4) -> bool:
    rcB = rc(seqB)
    tailA = seqA[-8:]
    tailB = rcB[-8:]
    matches = sum(1 for a, b in zip(tailA, tailB) if a == b)
    return matches >= min_match

def compute_primer_features(seq: str,
                            template: str,
                            all_primers: dict) -> dict:
    tm = p3.calc_tm(seq)
    hairpin_dG = compute_hairpin(seq)
    homodimer_dG = compute_homodimer(seq)
    mismatches = best_match_mismatches(seq, template)

    hetero = {}
    dangerous_3prime = {}
    for nameB, seqB in all_primers.items():
        if seqB == seq:
            continue
        dG = compute_heterodimer(seq, seqB)
        hetero[nameB] = dG
        dangerous_3prime[nameB] = has_3prime_dimer(seq, seqB)

    return {
        "tm": tm,
        "hairpin_dG": hairpin_dG,
        "homodimer_dG": homodimer_dG,
        "heterodimer_dG": hetero,
        "dangerous_3prime": dangerous_3prime,
        "mismatches": mismatches
    }

############################################################
# BIOPHYSICAL LAMP PRIMER SCORING
############################################################

def fraction_sequestered(dG_struct: float, scale: float = 4.0) -> float:
    dG = dG_struct / 100.0
    return 1.0 / (1.0 + math.exp(dG / scale))

def fraction_free(dG_struct: float, scale: float = 4.0) -> float:
    return 1.0 - fraction_sequestered(dG_struct, scale=scale)

def binding_efficiency(mismatches: int, per_mismatch_drop: float = 0.25) -> float:
    return max(0.0, (1.0 - per_mismatch_drop) ** mismatches)

def dimer_sink_factor(dangerous_3prime_flags: dict, per_dimer_drop: float = 0.5) -> float:
    n = sum(1 for v in dangerous_3prime_flags.values() if v)
    return (1.0 - per_dimer_drop) ** n

def lamp_tm_score(tm: float) -> float:
    return max(0.0, 1.0 - abs(tm - 60.0) / 12.0)

def lamp_primer_biophysical_score(feat: dict) -> float:
    f_free_struct = (
        fraction_free(feat["hairpin_dG"]) *
        fraction_free(feat["homodimer_dG"])
    )
    f_bind = binding_efficiency(feat["mismatches"])
    f_dimer = dimer_sink_factor(feat["dangerous_3prime"])
    tm_term = lamp_tm_score(feat["tm"])
    score = tm_term * f_free_struct * f_bind * f_dimer
    return max(0.0, score)

############################################################
# RT-LAMP SIMULATION
############################################################

def rt_lamp_rate(FIP_feat: dict,
                 BIP_feat: dict,
                 outer_feats: list) -> float:
    inner = lamp_primer_biophysical_score(FIP_feat) + lamp_primer_biophysical_score(BIP_feat)
    outer = sum(lamp_primer_biophysical_score(f) for f in outer_feats)
    combined = 0.7 * inner + 0.3 * outer
    base_k = 0.15
    scale_k = 0.40
    return base_k + scale_k * combined

def simulate_rt_qlamp(FIP_seq: str,
                      BIP_seq: str,
                      outer_seqs: list,
                      template_seq: str,
                      N0: float = 800,
                      t_max: float = 60,
                      dt: float = 0.25):

    all_primers = {
        "FIP": FIP_seq,
        "BIP": BIP_seq,
        "F3": outer_seqs[0],
        "F2": outer_seqs[1],
        "B3": outer_seqs[2]
    }

    FIP_feat = compute_primer_features(FIP_seq, template_seq, all_primers)
    BIP_feat = compute_primer_features(BIP_seq, template_seq, all_primers)
    outer_feats = [
        compute_primer_features(p, template_seq, all_primers)
        for p in outer_seqs
    ]

    k_eff = rt_lamp_rate(FIP_feat, BIP_feat, outer_feats)

    threshold_copies = 1e6
    plateau = 9e8
    t_rt = 8.0
    rt_eff = 0.45

    times, copies = [], []
    tt = None
    cDNA = 0.0
    t = 0.0

    while t <= t_max:
        if t < t_rt:
            cDNA = N0 * (t / t_rt) * rt_eff
            N = cDNA
        else:
            N = cDNA * math.exp(k_eff * (t - t_rt))

        if N > plateau:
            N = plateau

        times.append(t)
        copies.append(N)

        if tt is None and N >= threshold_copies:
            tt = t

        t += dt

    return np.array(times), np.array(copies), tt


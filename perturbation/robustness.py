"""
Robustness scoring for perturbation analysis.

Includes:
- robustness matrix construction
- perturbation robustness scoring
"""

import numpy as np
import primer3
from primers.mismatch import mismatch_cost
from thermodynamics.vienna_eval import vienna_fold_dna
from perturbation.mutate import perturb_primer


def build_robustness_matrix(lamp_sets, aln, names):
    roles = ["F3","F2","F1","B1","B2","B3"]
    matrix = []
    col_labels = []

    for idx, s in enumerate(lamp_sets):
        row_scores = []

        for r in roles:
            p = s[r]
            pos = s[f"{r}_pos"]
            variants = perturb_primer(p)

            for var in variants:
                matched = 0
                for i, seq in enumerate(aln):
                    seg = seq[pos:pos+len(var)]
                    if len(seg) != len(var):
                        continue
                    cost = mismatch_cost(var, seg)
                    if cost <= 4.0:
                        matched += 1

                ta = primer3.thermoanalysis.ThermoAnalysis()
                tm = ta.calc_tm(var)
                hp = ta.calc_hairpin(var)
                hd = ta.calc_homodimer(var)
                hairpin_tm = hp.tm if hp.structure_found else 0.0
                dimer_tm = hd.tm if hd.structure_found else 0.0

                v = vienna_fold_dna(var)
                mfe = v["mfe"]

                R = (1.0 * matched) + (0.5 * (tm - hairpin_tm - dimer_tm)) + (0.5 * (-mfe))

                row_scores.append(R)
                col_labels.append(f"{r}_{var}")

        matrix.append(row_scores)

    max_len = max(len(row) for row in matrix)
    mat_fixed = np.array([row + [np.nan]*(max_len - len(row)) for row in matrix])

    return mat_fixed, col_labels


def compute_perturbation_robustness_scores(mat, scored, threshold=10.0):
    scores = []

    for i in range(mat.shape[0]):
        row = mat[i, :]
        row = row[~np.isnan(row)]

        if len(row) == 0:
            scores.append({
                "index": i,
                "mean": np.nan,
                "min": np.nan,
                "frac_above": np.nan,
                "robust_score": np.nan
            })
            continue

        mean_val = np.mean(row)
        min_val  = np.min(row)
        frac_above = np.sum(row > threshold) / len(row)

        robust_score = (1.0 * mean_val) + (0.5 * min_val) + (0.5 * frac_above)

        scores.append({
            "index": i,
            "mean": mean_val,
            "min": min_val,
            "frac_above": frac_above,
            "robust_score": robust_score
        })

    for s, sc in zip(scored, scores):
        s["perturb_mean"]   = sc["mean"]
        s["perturb_min"]    = sc["min"]
        s["perturb_frac"]   = sc["frac_above"]
        s["perturb_score"]  = sc["robust_score"]

    return scores

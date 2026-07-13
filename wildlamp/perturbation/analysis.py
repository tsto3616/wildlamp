"""
Perturbation analysis for LAMP primer sets.

Evaluates:
- species coverage under mutation
- Primer3 thermodynamics
- ViennaRNA folding
"""

import primer3
from wildlamp.primers.mismatch import mismatch_cost
from wildlamp.thermodynamics.vienna_eval import vienna_fold_dna
from wildlamp.perturbation.mutate import perturb_primer


def perturbation_analysis(lset, aln, names, cost_max=4.0):
    roles = ["F3","F2","F1","B1","B2","B3"]
    results = []

    for r in roles:
        p = lset[r]
        pos = lset[f"{r}_pos"]

        variants = perturb_primer(p)

        for var in variants:
            matched = []
            failed = []

            # Species coverage
            for i, seq in enumerate(aln):
                seg = seq[pos:pos+len(var)]
                if len(seg) != len(var):
                    failed.append(names[i])
                    continue
                cost = mismatch_cost(var, seg)
                if cost <= cost_max:
                    matched.append(names[i])
                else:
                    failed.append(names[i])

            # Thermodynamics
            ta = primer3.thermoanalysis.ThermoAnalysis()
            tm = ta.calc_tm(var)
            hp = ta.calc_hairpin(var)
            hd = ta.calc_homodimer(var)

            # ViennaRNA
            v = vienna_fold_dna(var)

            results.append({
                "role": r,
                "original_primer": p,
                "perturbed_primer": var,
                "primer_pos": pos,
                "species_matched": len(matched),
                "species_failed": ",".join(failed),
                "tm": tm,
                "hairpin_tm": hp.tm if hp.structure_found else 0.0,
                "dimer_tm": hd.tm if hd.structure_found else 0.0,
                "mfe": v["mfe"],
                "structure": v["structure"]
            })

    return results

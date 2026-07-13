"""
Single‑base mutation utilities for perturbation analysis.
"""

def mutate_base(b: str):
    """Return all possible single-base mutations of a nucleotide."""
    return [x for x in "ACGT" if x != b]


def perturb_primer(primer: str):
    """Generate all single-base perturbations of a primer."""
    variants = []
    for i, b in enumerate(primer):
        for mb in mutate_base(b):
            variants.append(primer[:i] + mb + primer[i+1:])
    return variants

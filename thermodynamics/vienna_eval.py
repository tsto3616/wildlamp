"""
ViennaRNA folding and primer validation utilities for the wildlamp package.

Provides:
- vienna_fold_dna(): fold DNA (converted to RNA) using ViennaRNA Python API
- validate_primer_vienna(): evaluate individual primer stability
- validate_lamp_set_vienna(): evaluate all primers in a LAMP set
"""

import RNA


def vienna_fold_dna(dna_seq: str) -> dict:
    """
    Fold a DNA sequence using ViennaRNA (converted to RNA internally).

    Parameters
    ----------
    dna_seq : str
        DNA primer sequence.

    Returns
    -------
    dict
        {
            "rna_seq": str,
            "structure": str,
            "mfe": float
        }
    """
    rna_seq = dna_seq.replace("T", "U")
    structure, mfe = RNA.fold(rna_seq)
    return {
        "rna_seq": rna_seq,
        "structure": structure,
        "mfe": mfe
    }


def validate_primer_vienna(dna_seq: str, cutoff: float = -5.0) -> dict:
    """
    Validate a primer using ViennaRNA MFE folding.

    Parameters
    ----------
    dna_seq : str
        Primer sequence.
    cutoff : float
        Minimum acceptable MFE (higher = less stable structure).

    Returns
    -------
    dict
        {
            "primer": str,
            "rna_seq": str,
            "structure": str,
            "mfe": float,
            "pass": bool
        }
    """
    res = vienna_fold_dna(dna_seq)
    return {
        "primer": dna_seq,
        "rna_seq": res["rna_seq"],
        "structure": res["structure"],
        "mfe": res["mfe"],
        "pass": res["mfe"] > cutoff
    }


def validate_lamp_set_vienna(s: dict, cutoff: float = -5.0) -> dict:
    """
    Validate all primers in a LAMP set using ViennaRNA.

    Parameters
    ----------
    s : dict
        LAMP set containing keys:
        "F3","F2","F1","B3","B2","B1".
    cutoff : float
        Minimum acceptable MFE.

    Returns
    -------
    dict
        Modified LAMP set with:
        - s["vienna_validation"]: dict of per‑primer results
        - s["vienna_pass"]: bool
    """
    roles = ["F3", "F2", "F1", "B3", "B2", "B1"]
    results = {}
    all_pass = True

    for r in roles:
        res = validate_primer_vienna(s[r], cutoff=cutoff)
        results[r] = res
        if not res["pass"]:
            all_pass = False

    s["vienna_validation"] = results
    s["vienna_pass"] = all_pass
    return s

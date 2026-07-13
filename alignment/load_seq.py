from Bio import SeqIO

"""
Alignment utilities for wildlamp.

This module provides:
- load_alignment(): load FASTA sequences
- clean_alignment(): remove gap-only columns from multi-species alignments

These functions are used by the primer design pipeline.
"""

def load_alignment(path: str):
    """
    Load a FASTA alignment and return:
    - names: list of sequence IDs
    - seqs: list of uppercase DNA sequences (U→T)

    Parameters
    ----------
    path : str
        Path to FASTA alignment file.

    Returns
    -------
    names : list[str]
    seqs : list[str]
    """
    names, seqs = [], []
    for rec in SeqIO.parse(path, "fasta"):
        names.append(rec.id)
        seqs.append(str(rec.seq).upper().replace("U", "T"))
    return names, seqs


def clean_alignment(aln: list[str]):
    """
    Remove columns where all species have gaps ('-').

    Parameters
    ----------
    aln : list[str]
        List of aligned sequences.

    Returns
    -------
    cleaned : list[str]
        Alignment with gap-only columns removed.
    """
    max_len = max(len(s) for s in aln)
    padded = [s.ljust(max_len, "-") for s in aln]

    cols = list(zip(*padded))
    clean_cols = [col for col in cols if any(b in "ACGT" for b in col)]

    cleaned = ["".join(col[i] for col in clean_cols) for i in range(len(aln))]
    return cleaned

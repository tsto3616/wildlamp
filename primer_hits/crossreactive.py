#!/usr/bin/env python3

import csv
from Bio.Blast import NCBIWWW, NCBIXML

# ------------------------------------------------------------
# BLAST function
# ------------------------------------------------------------

def blast_sequence(seq, taxid):
    """
    Run BLASTN against nt restricted to a user-specified taxid.
    """
    result_handle = NCBIWWW.qblast(
        program="blastn",
        database="nt",
        sequence=seq,
        entrez_query=f"txid{taxid}[Organism]",
        expect=1000,
        word_size=7
    )
    return NCBIXML.read(result_handle)

# ------------------------------------------------------------
# Filtering logic
# ------------------------------------------------------------

def biologically_meaningful(hit, primer_len):
    """
    Apply biological filters: alignment length, mismatches, identity.
    """
    if hit.align_length < int(primer_len * 0.8):
        return False
    mismatches = hit.align_length - hit.identities
    if mismatches > 3:
        return False
    identity_pct = (hit.identities / hit.align_length) * 100
    if identity_pct < 85:
        return False
    return True

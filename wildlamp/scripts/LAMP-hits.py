#!/usr/bin/env python3

import csv
from Bio.Blast import NCBIWWW, NCBIXML

# ------------------------------------------------------------
# Panel definition (insert your real sequences)
# ------------------------------------------------------------

panel = {
    "DIO1": {
        "PrimerSet1": {
            "FIP": "TGGGAGGAGAGCCAAGACGTAGTCTGGAAAGATTCC",
            "BIP": "GCTAACTTGCTGCTTTTTTAGCTCAGATTGGTCTGCT",
        },
        "PrimerSet2": {
            "FIP": "TGGGAGGAGAGCCAAGACCTAGTTCTAGGTGATCAA",
            "BIP": "GCTAACTTGCTGCTTTTTTAGCTCAGATTGGTCTGCT"
        }
    },
    "DIO2": {
        "PrimerSet1": {
            "FIP": "TGCTGCACATCGGTCTTCCTTCTTTGTCTTTTGAAGTG",
            "BIP": "TCTTCTCCAGCCAACGCCGAAAGGGCCCCTTCTGCT"
        },
        "PrimerSet2": {
            "FIP": "TGCTGCACATCGGTCTTCCTTCTTTGTCTTTTGAAGTG",
            "BIP": "TCTTCTCCAGCCAACGCCGAAAGGGCCCCTTCTGCT"
        }
    },
    "SLC25A43": {
        "PrimerSet1": {
            "FIP": "GTAGAAAAAGCATGGAGGTGTGCAGAACCTGCTGGA",
            "BIP": "GAGAGCGTCTGGGTTACGCCTCCAGAACTTCGCCAA"
        },
        "PrimerSet2": {
            "FIP": "GTAGAAAAAGCATGGAGGTGTGCAGAACCTGCTGGA",
            "BIP": "GAGAGCGTCTGGGTTACGCCTCCAGAACTTCGCCAA"
        },    
    }, 
    "HSP70": {
        "PrimerSet1": {
            "FIP": "AGTGTCTGCAAAGGTCTTATGTCACAATAAAGCCATC",
            "BIP": "CCAAGTGTGACTACTGTAAAATGATAGTCCCCATATTAGTTTTGAATG"
        },
        "PrimerSet2": {
            "FIP": "AGTGTCTGCAAAGGTCTTATGTCACAATAAAGCCATC",
            "BIP": "ACTCGCAACACATACTAACCCCATATTAGTTTTGAATG"
        }
    },
}

# ------------------------------------------------------------
# BLAST function
# ------------------------------------------------------------

def blast_sequence(seq):
    result_handle = NCBIWWW.qblast(
        program="blastn",
        database="nt",
        sequence=seq,
        entrez_query=f"txid34884[Organism]",
        expect=1000,
        word_size=7
    )
    return NCBIXML.read(result_handle)

# ------------------------------------------------------------
# Filtering logic
# ------------------------------------------------------------

def biologically_meaningful(hit, primer_len):
    if hit.align_length < int(primer_len * 0.8):
        return False
    mismatches = hit.align_length - hit.identities
    if mismatches > 3:
        return False
    identity_pct = (hit.identities / hit.align_length) * 100
    if identity_pct < 85:
        return False
    return True

# ------------------------------------------------------------
# CSV output
# ------------------------------------------------------------

with open("callorhinus_crossreactivity.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow([
        "gene", "primer_set", "primer_name", "species",
        "scaffold", "start", "end",
        "alignment_length", "percent_identity", "mismatches", "pass_fail"
    ])

    for gene, sets in panel.items():
        for primer_set, primers in sets.items():
            for primer_name, seq in primers.items():

                blast_record = blast_sequence(seq)

                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:

                        primer_len = len(seq)
                        identity_pct = (hsp.identities / hsp.align_length) * 100
                        mismatches = hsp.align_length - hsp.identities

                        passed = biologically_meaningful(hsp, primer_len)

                        writer.writerow([
                            gene,
                            primer_set,
                            primer_name,
                            "Callorhinus ursinus",
                            alignment.title,
                            hsp.sbjct_start,
                            hsp.sbjct_end,
                            hsp.align_length,
                            f"{identity_pct:.2f}",
                            mismatches,
                            "PASS" if passed else "FAIL"
                        ])

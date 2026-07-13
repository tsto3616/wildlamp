"""
This command generates the primer hits against a genome using the taxid
"""

from primer_hits.crossreactive import blast_sequence, biologically_meaningful
import csv

def primers2hits(panel, taxid, output="crossreactivity_output.csv"):
    """
    Run BLAST crossreactivity on a list of primer sets.

    panel format:
        [
            {"FIP": "ATCG...", "BIP": "GGTA..."},
            {"FIP": "...", "BIP": "..."},
            ...
        ]

    Returns:
        output CSV filename
    """

    with open(output, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "primer_set_index", "primer_name", "sequence",
            "taxid", "scaffold", "start", "end",
            "alignment_length", "percent_identity", "mismatches", "pass_fail"
        ])

        for idx, primers in enumerate(panel):
            for primer_name, seq in primers.items():

                print(f"BLASTing set {idx} primer {primer_name} against taxid {taxid}...")

                blast_record = blast_sequence(seq, taxid)

                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:

                        primer_len = len(seq)
                        identity_pct = (hsp.identities / hsp.align_length) * 100
                        mismatches = hsp.align_length - hsp.identities

                        passed = biologically_meaningful(hsp, primer_len)

                        writer.writerow([
                            idx,
                            primer_name,
                            seq,
                            taxid,
                            alignment.title,
                            hsp.sbjct_start,
                            hsp.sbjct_end,
                            hsp.align_length,
                            f"{identity_pct:.2f}",
                            mismatches,
                            "PASS" if passed else "FAIL"
                        ])

    print(f"Crossreactivity table written to {output}")
    return output

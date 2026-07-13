"""
This is for the checking of primers from a csv file to the genome specified by the user through a 
ncbi taxid search
"""

from wildlamp.pipeline.primer2hits import primers2hits

def csv2hits(csv_file, taxid, output="crossreactivity_output.csv"):
    """
    Read primer sets from a CSV and run BLAST.

    CSV format:
        fip,bip
        TGGGAGGAGAGCCAAGACGTAGTCTGGAAAGATTCC,GCTAACTTGCTGCTTTTTTAGCTCAGATTGGTCTGCT
        TGCTGCACATCGGTCTTCCTTCTTTGTCTTTTGAAGTG,TCTTCTCCAGCCAACGCCGAAAGGGCCCCTTCTGCT

    Returns:
        output CSV filename
    """

    panel = []

    with open(csv_file, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            fip = row["fip"].strip()
            bip = row["bip"].strip()
            panel.append({"FIP": fip, "BIP": bip})

    return primers2hits(panel, taxid=taxid, output=output)

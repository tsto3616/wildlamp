import csv

def load_gene_primer_csv(csv_path):
    """
    CSV must contain columns:
        Gene, fasta, F1, F2, B1, B2, F3, B3
    """
    genes = {}
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            genes[row["Gene"]] = {
                "fasta": row["fasta"],
                "F1": row["F1"],
                "F2": row["F2"],
                "B1": row["B1"],
                "B2": row["B2"],
                "F3": row["F3"],
                "B3": row["B3"]
            }
    return genes


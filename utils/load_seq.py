import csv
from Bio import SeqIO

def load_primers_csv(path):
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows

def load_fasta(path):
    rec = next(SeqIO.parse(path, "fasta"))
    return str(rec.seq).upper()

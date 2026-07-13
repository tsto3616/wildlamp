import csv

def load_primers_from_csv(csv_path):
    primers = {}
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            role = row["role"].strip()
            seq = row["sequence"].strip()
            primers[role] = seq
    return primers

"""
Accession extraction utilities for wildlamp.

Provides:
- extract_accession(): parse accession IDs from text
- extract_accessions_from_csv(): read species_list column and extract accessions
"""

import re
import csv


def extract_accession(text: str) -> str:
    """Extract RefSeq-style or generic accession from text."""
    m = re.search(r"(X[MP]_\d+\.\d+)", text)
    if m:
        return m.group(1)

    m = re.search(r"([A-Z]{2}_\d+)", text)
    if m:
        return m.group(1)

    m = re.search(r"ref\|(.*?)\|", text)
    if m:
        return m.group(1)

    return text.strip()


def extract_accessions_from_csv(csv_file: str):
    """Extract all accession IDs from species_list column of a CSV."""
    accessions = set()

    with open(csv_file, newline="") as f:
        reader = csv.DictReader(f)

        for row in reader:
            if "species_list" not in row:
                continue

            species_entries = row["species_list"].split(",")

            for entry in species_entries:
                acc = extract_accession(entry)
                accessions.add(acc)

    return sorted(accessions)

"""
Species lookup table writer for wildlamp.

Provides:
- write_species_lookup(): write accession → species → gene symbol table
"""

import csv
from wildlamp.species.ncbi_lookup import get_species_from_accession


def write_species_lookup(accessions, outfile="species_lookup.csv"):
    """Write a CSV table mapping accession → species → gene symbol."""
    with open(outfile, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Accession", "Species_name", "Gene_symbol"])

        for acc in accessions:
            species, gene = get_species_from_accession(acc)
            writer.writerow([
                acc,
                species if species else "UNKNOWN",
                gene if gene else "UNKNOWN"
            ])

    print(f"Species lookup written to {outfile}")

"""
Pipeline wrapper for species lookup from validated primer CSV.
"""

from ..species.accession import extract_accessions_from_csv
from ..species.table import write_species_lookup


def species_fetch(csv_file):
    """Extract accessions from a validated primer CSV and write species lookup."""
    accessions = extract_accessions_from_csv(csv_file)
    write_species_lookup(accessions)
    print("Species_fetch finished")

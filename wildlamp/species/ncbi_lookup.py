"""
NCBI lookup utilities for wildlamp.

Provides:
- accession_to_gene_id(): nuccore → gene ID
- get_species_from_accession(): gene ID → species + gene symbol
"""

import requests


def accession_to_gene_id(acc: str):
    """Convert a nuccore accession to a Gene ID using NCBI elink."""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
    params = {
        "dbfrom": "nuccore",
        "db": "gene",
        "id": acc,
        "retmode": "json"
    }

    try:
        r = requests.get(url, params=params, timeout=10)
        data = r.json()

        linksets = data.get("linksets", [])
        if not linksets:
            return None

        links = linksets[0].get("linksetdbs", [])
        if not links:
            return None

        gene_ids = links[0].get("links", [])
        if not gene_ids:
            return None

        return gene_ids[0]

    except Exception:
        return None


def get_species_from_accession(acc: str):
    """Return (species_name, gene_symbol) for a nuccore accession."""
    gene_id = accession_to_gene_id(acc)
    if gene_id is None:
        return None, None

    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    summary_params = {
        "db": "gene",
        "id": gene_id,
        "retmode": "json"
    }

    try:
        r = requests.get(summary_url, params=summary_params, timeout=10)
        data = r.json()

        result = data["result"][str(gene_id)]

        species_name = result.get("organism")
        gene_symbol = result.get("name")

        return species_name, gene_symbol

    except Exception:
        return None, None

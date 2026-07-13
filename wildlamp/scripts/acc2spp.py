import csv
import re
import requests

# this script was used to label the species for the accessions and generate a species table foe each primer

def extract_accession(text):
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


def accession_to_gene_id(acc):
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

# query against NCBI to get species name and gene symbol for a given accession
def get_species_from_accession(acc):

    gene_id = accession_to_gene_id(acc)
    if gene_id is None:
        return None, None

    # Fetch gene summary
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


# cross check the accessions against the CSV file and extract accessions
def extract_accessions_from_csv(csv_file):
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


# write as a csv table
def write_species_lookup(accessions, outfile="species_lookup.csv"):
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


# run it
if __name__ == "__main__":
    csv_file = "DOI3_rtlamp_validated_top20.csv"

    accessions = extract_accessions_from_csv(csv_file)
    write_species_lookup(accessions)

from wildlamp.utils.load_seq import load_primers_csv, load_fasta
from wildlamp.thermodynamics.fold_rna import fold_rna
from wildlamp.primers.analyse_primer_rows import analyse_row
from wildlamp.utils.valid_csv_write import write_outputs

def final_valid(primer_csv, fasta_path, output="fill with csv"):
    """
    CENTRAL USER COMMAND.
    - Loads primers from CSV (each row independent)
    - Loads FASTA template sequence
    - Computes mismatch, accessibility, ΔG, GC%, MFE
    - Writes unified CSV output
    """

    template_seq = load_fasta(fasta_path)
    template_structure, template_mfe = fold_rna(template_seq)

    primer_rows = load_primers_csv(primer_csv)

    results = []
    for row in primer_rows:
        result = analyse_row(row, template_seq, template_structure, template_mfe)
        results.append(result)

    write_outputs(results, output)

    return results

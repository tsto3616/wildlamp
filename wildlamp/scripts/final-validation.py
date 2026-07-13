#!/usr/bin/env python

import RNA
import primer3.bindings as p3
from Bio import SeqIO
import csv

############################################################
# GENE CONFIGURATION
############################################################

GENES = {
    "DIO2": {
        "fasta": "DIO2_N_fur.fasta",
        "F3": "AGGGGCTGACTTCCTGTT",
        "F2": "AGGGGCTGACTTCCTGTT",
        "F1": "CAGGAAGACGGATGTGCA",
        "B3": "GCCTTTGAACCTGTGTGT",
        "B2": "GCCCCTTCTCCTACAACC",
        "B1": "ACGTTGGCTGGAGAAGAA"
    },
    "DIO1": {
        "fasta": "DIO1_N_fur.fasta",
        "F3": "TTCTGGAAAAGCTCCACA",
        "F2": "GTAGTCTGGAAAGATTCC",
        "F1": "CTTGACTCTCCTCCCAGC",
        "B3": "GACTCCCGTCCCCTGCTGA",
        "B2": "AGCTCGGATTGGTCTGCT",
        "B1": "AAAGCCGCAAGTTAGCTGA"
    },
    "DIO3": {
        "fasta": "DIO3_N_fur.fasta",
        "F3": "ACTCTCCCTCCCCTCCCT",
        "F2": "GCTGGCTCGACATCCCTT",
        "F1": "GAGTGTGGATGTGATGCC",
        "B3": "CCAGACGCAGAGAGACTT",
        "B2": "GCCTGCCGCTCGCCCTAGTT",
        "B1": "TTGTTCCTGAGTTGGGGGTG"
    },
    "SLC25A43": {
        "fasta": "slc25a43_N_fur.fasta",
        "F3": "GGCTGGGAGTCTTGCAGG",
        "F2": "TGTGCAGAACCTGCTGGA",
        "F1": "CCTCCATGCTTTTTCTAC",
        "B3": "GAAACTCTGGAATGGACG",
        "B2": "CCTCCAGAACTTCGCCAA",
        "B1": "CGTAACCCAGACGCTCTC"
    },
    "HSP70": {
        "fasta": "HSP70_N_fur.fasta",
        "F3": "TCTGCTTGGATTTCCGTG",
        "F2": "ATGTCACAATAAAGCCATC",
        "F1": "AAGACCTTTGCAGACACT",
        "B3": "GTAGCTTGCATTTGAACT",
        "B2": "CCCCATATTAGTTTTGAATG",
        "B1": "ACTATCATTTTACAGTAGTCACACTTGG"
    }
}

############################################################
# UTILITIES
############################################################

def load_gene(path):
    rec = next(SeqIO.parse(path, "fasta"))
    return str(rec.seq).upper()

def gc_content(seq):
    gc = seq.count("G") + seq.count("C")
    return (gc / len(seq)) * 100

def fold_rna(sequence):
    """Return dot-bracket structure using ViennaRNA."""
    fc = RNA.fold_compound(sequence.replace("T", "U"))
    structure, mfe = fc.mfe()
    return structure, mfe

def best_match(primer, template):
    """Return (mismatches, best_position)."""
    Lp = len(primer)
    Lt = len(template)
    best_mm = Lp
    best_pos = -1
    for i in range(Lt - Lp + 1):
        window = template[i:i+Lp]
        mm = sum(1 for a, b in zip(primer, window) if a != b)
        if mm < best_mm:
            best_mm = mm
            best_pos = i
            if mm == 0:
                break
    return best_mm, best_pos

def binding_accessibility(structure, start, end):
    region = structure[start:end]
    unpaired = region.count('.')
    return (unpaired / len(region)) * 100

def thermo(primer):
    hp = p3.calc_hairpin(primer).dg
    hd = p3.calc_homodimer(primer).dg
    return hp, hd

############################################################
# MAIN PER-GENE ANALYSIS
############################################################

def analyse_gene(gene_name, cfg): 
    gene_seq = load_gene(cfg["fasta"])
    rna_structure, mfe = fold_rna(gene_seq)

    primers = {
        "F3": cfg["F3"],
        "F2": cfg["F2"],
        "F1": cfg["F1"],
        "B3": cfg["B3"],
        "B2": cfg["B2"],
        "B1": cfg["B1"]
    }

    rows = []

    for pname, pseq in primers.items():
        mm, pos = best_match(pseq, gene_seq)

        if pos != -1:
            acc = binding_accessibility(rna_structure, pos, pos + len(pseq))
        else:
            acc = 0.0

        hp, hd = thermo(pseq)

        rows.append({
            "Gene": gene_name,
            "Primer": pname,
            "Sequence": pseq,
            "GC%": gc_content(pseq),
            "Mismatches": mm,
            "Binding_position": pos,
            "Accessibility_%": acc,
            "Hairpin_dG": hp,
            "Homodimer_dG": hd,
            "Transcript_MFE": mfe
        })

    return rows

############################################################
# WRITE TABLES
############################################################

def write_gene_tables(all_rows):
    with open("primer_binding_accessibility_thermo.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "Gene","Primer","Sequence","GC%","Mismatches",
            "Binding_position","Accessibility_%",
            "Hairpin_dG","Homodimer_dG","Transcript_MFE"
        ])
        w.writeheader()
        w.writerows(all_rows)

def write_gc_table(all_rows):
    with open("primer_gc_table.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["Gene","Primer","Sequence","GC%"])
        w.writeheader()
        for r in all_rows:
            w.writerow({
                "Gene": r["Gene"],
                "Primer": r["Primer"],
                "Sequence": r["Sequence"],
                "GC%": r["GC%"],
                
            })

############################################################
# WRITE LONG-FORM CSV FOR GC, Hairpin_dG, Homodimer_dG, Accessibility
############################################################

def write_long_all(all_rows, filename="primer_long_format.csv"):
    fieldnames = [
        "Gene",
        "Primer",
        "Sequence",
        "GC%",
        "Hairpin_dG",
        "Homodimer_dG",
        "Accessibility_%"
    ]

    with open(filename, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for row in all_rows:
            writer.writerow({
                "Gene": row["Gene"],
                "Primer": row["Primer"],
                "Sequence": row["Sequence"],
                "GC%": row["GC%"],
                "Hairpin_dG": row["Hairpin_dG"],
                "Homodimer_dG": row["Homodimer_dG"],
                "Accessibility_%": row["Accessibility_%"]
            })

    print(f"Long-format primer CSV written to: {filename}")

############################################################
# MAIN
############################################################

def main():
    all_rows = []

    for gene_name, cfg in GENES.items():
        print(f"\n=== ANALYSING {gene_name} ===")
        rows = analyse_gene(gene_name, cfg)
        all_rows.extend(rows)

        for r in rows:
            print(
                f"{r['Primer']:>3} | GC={r['GC%']:5.1f}% | "
                f"mm={r['Mismatches']:2d} | pos={r['Binding_position']:5d} | "
                f"acc={r['Accessibility_%']:5.1f}% | "
                f"hp_dG={r['Hairpin_dG']:7.2f} | hd_dG={r['Homodimer_dG']:7.2f}"
            )

    write_gene_tables(all_rows)
    write_gc_table(all_rows)
    write_long_all(all_rows)   # <-- NEW LINE

if __name__ == "__main__":
    main()

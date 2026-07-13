import os
import pandas as pd
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# USER INPUT: folder containing FASTA files
# ---------------------------------------------------------
fasta_folder = "fasta_folder/"   # change this once

# ---------------------------------------------------------
# INCLUSIVE MISMATCH SCORING
# ---------------------------------------------------------
def best_mismatch_score(primer, transcript):
    """
    Inclusive mismatch scoring:
    - Always returns a finite score
    - Allows mismatches anywhere
    - Penalises 3' mismatches strongly
    - Penalises internal mismatches moderately
    - Penalises 5' mismatches weakly
    """
    Lp = len(primer)
    Lt = len(transcript)

    if Lp == 0 or Lt < Lp:
        return 999  # worst-case but finite

    best_score = 999

    for i in range(Lt - Lp + 1):
        window = transcript[i:i+Lp]
        mismatches = 0
        penalty = 0

        for pos, (p, w) in enumerate(zip(primer, window)):
            if p != w:
                mismatches += 1
                if pos >= Lp - 5:      # 3' end
                    penalty += 3
                elif pos >= 5:         # internal
                    penalty += 2
                else:                  # 5' end
                    penalty += 1

        score = mismatches + penalty
        if score < best_score:
            best_score = score

    return best_score

# ---------------------------------------------------------
# AUTO-DETECT FASTA FILES
# ---------------------------------------------------------
fasta_files = [
    os.path.join(fasta_folder, f)
    for f in os.listdir(fasta_folder)
    if f.endswith(".fasta")
]

print(f"Found {len(fasta_files)} FASTA files:")
for f in fasta_files:
    print("  ", f)

# ---------------------------------------------------------
# CLASSIFY FILES: reference vs multi-transcript
# ---------------------------------------------------------
gene_fastas = {}

for f in fasta_files:
    records = list(SeqIO.parse(f, "fasta"))
    if len(records) == 0:
        continue

    gene_name = os.path.basename(f).split("_")[0]

    if gene_name not in gene_fastas:
        gene_fastas[gene_name] = {"ref": None, "multi": None}

    if len(records) == 1:
        gene_fastas[gene_name]["ref"] = f
    else:
        gene_fastas[gene_name]["multi"] = f

print("\nGene FASTA classification:")
for gene, files in gene_fastas.items():
    print(f"Gene {gene}: ref={files['ref']}, multi={files['multi']}")

# ---------------------------------------------------------
# LOAD PRIMER CSV
# ---------------------------------------------------------
df = pd.read_csv("primer_long_format.csv")

required_cols = ["Gene", "Primer", "Sequence", "Hairpin_dG", "Homodimer_dG", "Accessibility_%"]
for col in required_cols:
    if col not in df.columns:
        raise ValueError(f"Missing column in CSV: {col}")

print("\nLoaded primer CSV with columns:")
print(df.columns)

# ---------------------------------------------------------
# HEAVY PERTURBATION PARAMETERS
# ---------------------------------------------------------
N = 500
gc_drift = 0.03
dg_noise = 10
access_noise = 0.15

heavy_results = []

# ---------------------------------------------------------
# HEAVY PERTURBATION (REFERENCE TRANSCRIPT)
# ---------------------------------------------------------
for gene, files in gene_fastas.items():
    ref_file = files["ref"]
    if ref_file is None:
        print(f"\n[Heavy] Skipping gene {gene}: no reference FASTA.")
        continue

    ref_record = list(SeqIO.parse(ref_file, "fasta"))[0]
    ref_seq = str(ref_record.seq)

    gene_df = df[df["Gene"] == gene]
    if gene_df.empty:
        print(f"\n[Heavy] No primers for gene {gene}.")
        continue

    print(f"\n[Heavy] Processing gene {gene} on reference transcript {ref_record.id}.")

    for idx, row in gene_df.iterrows():
        primer = row["Primer"]
        primer_seq = str(row["Sequence"])

        mismatch_score = best_mismatch_score(primer_seq, ref_seq)

        base_gc = (primer_seq.count("G") + primer_seq.count("C")) / len(primer_seq)
        base_access = float(row["Accessibility_%"])
        base_hairpin = float(row["Hairpin_dG"])
        base_homodimer = float(row["Homodimer_dG"])

        deviations = []

        for _ in range(N):
            gc_p = base_gc + np.random.uniform(-gc_drift, gc_drift)
            hairpin_p = base_hairpin + np.random.uniform(-dg_noise, dg_noise)
            homodimer_p = base_homodimer + np.random.uniform(-dg_noise, dg_noise)
            access_p = base_access + np.random.uniform(-access_noise, access_noise)
            access_p = max(0, min(access_p, 1))

            gc_dev = abs(gc_p - base_gc) / base_gc if base_gc != 0 else abs(gc_p - base_gc)
            hairpin_dev = abs(hairpin_p - base_hairpin) / abs(base_hairpin) if base_hairpin != 0 else abs(hairpin_p - base_hairpin)
            homodimer_dev = abs(homodimer_p - base_homodimer) / abs(base_homodimer) if base_homodimer != 0 else abs(homodimer_p - base_homodimer)

            if base_access == 0:
                access_dev = abs(access_p - base_access)
            else:
                access_dev = abs(access_p - base_access) / base_access

            mismatch_p = mismatch_score + np.random.uniform(-0.5, 0.5)
            mismatch_p = max(0, mismatch_p)

            if mismatch_score == 0:
                mismatch_dev = mismatch_p
            else:
                mismatch_dev = abs(mismatch_p - mismatch_score) / mismatch_score

            global_dev = 100 * (gc_dev + hairpin_dev + homodimer_dev + access_dev + mismatch_dev) / 5
            deviations.append(global_dev)

        heavy_results.append({
            "Gene": gene,
            "Primer": primer,
            "MedianDeviation": np.median(deviations),
            "StdDeviation": np.std(deviations)
        })

heavy_df = pd.DataFrame(heavy_results)

print("\nHEAVY PERTURBATION RESULTS:")
print(heavy_df)

# ---------------------------------------------------------
# LIGHTWEIGHT MULTI-TRANSCRIPT BASELINE
# ---------------------------------------------------------
light_results = []

for gene, files in gene_fastas.items():
    multi_file = files["multi"]
    if multi_file is None:
        print(f"\n[Light] Skipping gene {gene}: no multi-transcript FASTA.")
        continue

    records = list(SeqIO.parse(multi_file, "fasta"))
    gene_df = df[df["Gene"] == gene]

    print(f"\n[Light] Processing gene {gene} on {len(records)} transcripts.")

    for idx, row in gene_df.iterrows():
        primer = row["Primer"]
        primer_seq = str(row["Sequence"])

        base_gc = (primer_seq.count("G") + primer_seq.count("C")) / len(primer_seq)
        base_access = float(row["Accessibility_%"])
        base_hairpin = float(row["Hairpin_dG"])
        base_homodimer = float(row["Homodimer_dG"])

        for rec in records:
            seq = str(rec.seq)
            mismatch_score = best_mismatch_score(primer_seq, seq)

            light_results.append({
                "Gene": gene,
                "Primer": primer,
                "Transcript": rec.id,
                "GC": base_gc,
                "Accessibility_%": base_access,
                "Hairpin_dG": base_hairpin,
                "Homodimer_dG": base_homodimer,
                "MismatchScore": mismatch_score
            })

light_df = pd.DataFrame(light_results)

print("\nLIGHTWEIGHT RESULTS:")
print(light_df)

light_df.to_csv("lightweight_table.csv")
# ---------------------------------------------------------
# PLOTS
# ---------------------------------------------------------

# HEAVY PLOT
if not heavy_df.empty:
    heavy_gene = (
        heavy_df.groupby("Gene")[["MedianDeviation", "StdDeviation"]]
        .mean()
        .reset_index()
    )

    plt.figure(figsize=(10, 6))
    plt.bar(
        heavy_gene["Gene"],
        heavy_gene["MedianDeviation"],
        yerr=heavy_gene["StdDeviation"],
        capsize=8,
        color="Red"
    )
    plt.ylabel("Perturbation Deviation (%)")
    plt.title("Heavy Monte‑Carlo Perturbation (Reference Transcript)")
    plt.grid(axis="y", linestyle="--", alpha=0.6)
    plt.tight_layout()
    plt.show()

# LIGHT PLOT (NO X-AXIS)
if not light_df.empty:
    plt.figure(figsize=(12, 6))

    for gene in df["Gene"].unique():
        subset = light_df[light_df["Gene"] == gene]
        if subset.empty:
            continue

        plt.scatter(
            range(len(subset)),
            subset["MismatchScore"],
            label=gene,
            s=60
        )

    plt.ylabel("Mismatch Score (lower = better)")
    plt.title("Lightweight Multi‑Transcript Mismatch Scan")

    plt.xticks([])
    plt.gca().set_xticks([])
    plt.gca().set_xlabel("")

    plt.legend()
    plt.grid(axis="y", linestyle="--", alpha=0.6)
    plt.tight_layout()
    plt.show()

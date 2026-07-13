import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from primers.primer_from_csv import load_primers_from_csv


# ---------------------------------------------------------
# Inclusive mismatch scoring
# ---------------------------------------------------------
def best_mismatch_score(primer, transcript):
    Lp = len(primer)
    Lt = len(transcript)

    if Lp == 0 or Lt < Lp:
        return 999

    best_score = 999

    for i in range(Lt - Lp + 1):
        window = transcript[i:i+Lp]
        mismatches = 0
        penalty = 0

        for pos, (p, w) in enumerate(zip(primer, window)):
            if p != w:
                mismatches += 1
                if pos >= Lp - 5:
                    penalty += 3
                elif pos >= 5:
                    penalty += 2
                else:
                    penalty += 1

        score = mismatches + penalty
        best_score = min(best_score, score)

    return best_score


# ---------------------------------------------------------
# Gene‑agnostic hidden command
# ---------------------------------------------------------
def find_best_matches(fasta_folder, primer_csv):
    """
    Gene‑agnostic version.
    - Reads primers using load_primers_from_csv()
    - No Gene column required
    - Every primer tested against every FASTA file
    - Heavy perturbation uses FIRST sequence in each FASTA
    - Light mismatch scan uses ALL sequences in each FASTA

    Returns:
        heavy_df, light_df
    """

    # ---------------------------------------------------------
    # Detect FASTA files
    # ---------------------------------------------------------
    fasta_files = [
        os.path.join(fasta_folder, f)
        for f in os.listdir(fasta_folder)
        if f.endswith(".fasta")
    ]

    if not fasta_files:
        raise ValueError("No FASTA files found in folder.")

    # ---------------------------------------------------------
    # Load primers from CSV (same format as final_primers)
    # ---------------------------------------------------------
    primers = load_primers_from_csv(primer_csv)

    required = ["F3","F2","F1","B3","B2","B1"]
    for r in required:
        if r not in primers:
            raise ValueError(f"Primer {r} missing from CSV.")

    # Convert to DataFrame for convenience
    primer_rows = []
    for role, seq in primers.items():
        primer_rows.append({
            "Primer": role,
            "Sequence": seq,
            "Hairpin_dG": 0.0,
            "Homodimer_dG": 0.0,
            "Accessibility_%": 1.0
        })
    df = pd.DataFrame(primer_rows)

    # ---------------------------------------------------------
    # Heavy perturbation parameters
    # ---------------------------------------------------------
    N = 500
    gc_drift = 0.03
    dg_noise = 10
    access_noise = 0.15

    heavy_results = []
    light_results = []

    # ---------------------------------------------------------
    # Loop over FASTA files
    # ---------------------------------------------------------
    for fasta_path in fasta_files:
        fasta_name = os.path.basename(fasta_path)

        records = list(SeqIO.parse(fasta_path, "fasta"))
        if not records:
            continue

        # Heavy perturbation uses FIRST sequence
        ref_seq = str(records[0].seq)

        # Light mismatch scan uses ALL sequences
        multi_seqs = [str(rec.seq) for rec in records]

        # ---------------------------------------------------------
        # Loop over primers
        # ---------------------------------------------------------
        for _, row in df.iterrows():
            primer = row["Primer"]
            primer_seq = str(row["Sequence"])

            base_gc = (primer_seq.count("G") + primer_seq.count("C")) / len(primer_seq)
            base_access = float(row["Accessibility_%"])
            base_hairpin = float(row["Hairpin_dG"])
            base_homodimer = float(row["Homodimer_dG"])

            # ---------------------------------------------------------
            # Heavy perturbation
            # ---------------------------------------------------------
            mismatch_score = best_mismatch_score(primer_seq, ref_seq)

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
                access_dev = abs(access_p - base_access) / base_access if base_access != 0 else abs(access_p - base_access)

                mismatch_p = mismatch_score + np.random.uniform(-0.5, 0.5)
                mismatch_p = max(0, mismatch_p)

                mismatch_dev = mismatch_p if mismatch_score == 0 else abs(mismatch_p - mismatch_score) / mismatch_score

                global_dev = 100 * (gc_dev + hairpin_dev + homodimer_dev + access_dev + mismatch_dev) / 5
                deviations.append(global_dev)

            heavy_results.append({
                "FASTA": fasta_name,
                "Primer": primer,
                "MedianDeviation": np.median(deviations),
                "StdDeviation": np.std(deviations)
            })

            # ---------------------------------------------------------
            # Light mismatch scan
            # ---------------------------------------------------------
            for rec in records:
                seq = str(rec.seq)
                mismatch_score = best_mismatch_score(primer_seq, seq)

                light_results.append({
                    "FASTA": fasta_name,
                    "Primer": primer,
                    "Transcript": rec.id,
                    "GC": base_gc,
                    "Accessibility_%": base_access,
                    "Hairpin_dG": base_hairpin,
                    "Homodimer_dG": base_homodimer,
                    "MismatchScore": mismatch_score
                })

    heavy_df = pd.DataFrame(heavy_results)
    light_df = pd.DataFrame(light_results)

    return heavy_df, light_df

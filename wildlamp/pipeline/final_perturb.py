# wildlamp/perturbation/final_perturbation.py

import pandas as pd
import matplotlib.pyplot as plt
from wildlamp.perturbation.best_mismatch_score import find_best_matches


def final_perturbation(fasta_folder, primer_csv, output="perturbation_output.csv"):
    """
    Public command.
    Runs:
    - heavy perturbation
    - lightweight mismatch scan
    - writes CSV
    - plots results

    Gene‑agnostic version:
    - No Gene column required
    - Grouping is done by FASTA file and Primer name
    """

    heavy_df, light_df = find_best_matches(fasta_folder, primer_csv)

    # ---------------------------------------------------------
    # Write raw tables
    # ---------------------------------------------------------
    heavy_df.to_csv("heavy_table.csv", index=False)
    light_df.to_csv("light_table.csv", index=False)

    # ---------------------------------------------------------
    # Combine heavy + light (gene‑agnostic)
    # ---------------------------------------------------------
    # Light mismatch scores averaged per FASTA + Primer
    light_summary = (
        light_df.groupby(["FASTA", "Primer"])["MismatchScore"]
        .mean()
        .reset_index()
        .rename(columns={"MismatchScore": "MeanMismatchScore"})
    )

    combined = heavy_df.merge(
        light_summary,
        on=["FASTA", "Primer"],
        how="left"
    )

    combined.to_csv(output, index=False)
    print(f"\nSaved combined perturbation table to {output}")

    # ---------------------------------------------------------
    # Heavy plot (per FASTA)
    # ---------------------------------------------------------
    if not heavy_df.empty:
        heavy_fasta = (
            heavy_df.groupby("FASTA")[["MedianDeviation", "StdDeviation"]]
            .mean()
            .reset_index()
        )

        plt.figure(figsize=(10, 6))
        plt.bar(
            heavy_fasta["FASTA"],
            heavy_fasta["MedianDeviation"],
            yerr=heavy_fasta["StdDeviation"],
            capsize=8,
            color="Red"
        )
        plt.ylabel("Perturbation Deviation (%)")
        plt.title("Heavy Monte‑Carlo Perturbation (Reference Transcript)")
        plt.grid(axis="y", linestyle="--", alpha=0.6)
        plt.tight_layout()
        plt.show()

    # ---------------------------------------------------------
    # Light plot (per FASTA)
    # ---------------------------------------------------------
    if not light_df.empty:
        plt.figure(figsize=(12, 6))

        for fasta_name in light_df["FASTA"].unique():
            subset = light_df[light_df["FASTA"] == fasta_name]
            plt.scatter(
                range(len(subset)),
                subset["MismatchScore"],
                label=fasta_name,
                s=60
            )

        plt.ylabel("Mismatch Score (lower = better)")
        plt.title("Lightweight Multi‑Transcript Mismatch Scan")
        plt.xticks([])
        plt.legend()
        plt.grid(axis="y", linestyle="--", alpha=0.6)
        plt.tight_layout()
        plt.show()

    return {
        "heavy": heavy_df,
        "light": light_df,
        "output_csv": output
    }
